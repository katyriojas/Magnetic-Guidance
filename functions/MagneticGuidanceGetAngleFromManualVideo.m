%% MagneticGuidanceGetAngleFromVideo
%
%   insertion_angle.time  => [s] time since start of insertion
%   insertion_angle.angle => [deg] angle of the tip in cochlea frame
%   insertion_angle.angle_smooth => angle smoothed by smooth_span
%
% Trevor Bruns and Katy Riojas
% September 2019

% Sample video location
clear all; close all; clc;
folderpath_vid = 'C:\Users\riojaske\Documents\magsteer\Magnetic-Guidance\data\phantom\manual';

% filename_vid = '\videos\trial1\trial1-manual-ea1-cropped-60fps.mp4';
% filename_save = strcat(folderpath_vid,'\pman1_60fps_angular_depth.mat');
% filename_vid = '\videos\trial2\trial2-manual-ea1-cropped-60fps.mp4';
% filename_save = strcat(folderpath_vid,'\pman2_60fps_angular_depth.mat');
% filename_vid = '\videos\trial3\trial3-manual-ea1-cropped-60fps.mp4';
% filename_save = strcat(folderpath_vid,'\pman3_60fps_angular_depth.mat');
filename_vid = '\videos\trial4\trial4-manual-ea2-cropped-60fps_divide.mp4';
filename_save = strcat(folderpath_vid,'\pman4_60fps_angular_depth.mat');

video_path = fullfile(folderpath_vid, filename_vid);

% This should be equal to robotic phantom data if the frame rates are the same
smooth_span = 40; % number of samples to smooth
magnification = 300;

% Load video and read first frame
vid = VideoReader(video_path);
start_time = vid.CurrentTime;
curr_frame = readFrame(vid);
imshow(curr_frame,'InitialMagnification',magnification); % show frame

% Center Point
fprintf('Select Center Location\n');
[centerX,centerY] = ginput(1);
centerPnt = [centerX;centerY];

% 0 degree marker
fprintf('Select 0 degree Mark\n');
[angle0X,angle0Y] = ginput(1);
angle0Pnt = [angle0X;angle0Y];
clf();

% 0 degree vector
angle0.vec = angle0Pnt - centerPnt;

%% Step through video and determine angle
num_frames = floor(vid.FrameRate*vid.Duration); % Total frames
insertion_angle.time  = zeros(1,num_frames); % [s]
insertion_angle.angle = zeros(1,num_frames); % [deg]
frame_count = 0;
vid.CurrentTime = 0;

while hasFrame(vid)
    if frame_count == 0
      fprintf('Select Tip location (most distal end of most apical electrode pad\n');
    end
    
    frame_count = frame_count + 1;
    insertion_angle.time(frame_count) = vid.CurrentTime;
    curr_frame = readFrame(vid); % read next frame, this will increment time
    imshow(curr_frame,'InitialMagnification',magnification); % show frame
    
    percDone = 100*frame_count/num_frames;
    fprintf('%2.3f perc\n',percDone)
    
    % Segment tip marker
    [tipX,tipY] = ginput(1);
    tipPnt = [tipX;tipY];
    
    % Find angle
    tip.vec = tipPnt - centerPnt;
    insertion_angle.angle(frame_count) = rad2deg(vectorAngle(angle0.vec',tip.vec'));
    
end

% Ensure time increases monotonically (last frame shows previous time for some reason)
%delete after the fact if we want to.
% insertion_angle.time(end)  = [];
% insertion_angle.angle(end) = [];

%% Account 360 degree 'rollover'
ang_temp = insertion_angle.angle;
for ii = 2:length(ang_temp)
    if (ang_temp(ii)-ang_temp(ii-1)) > 180
        ang_temp = ang_temp - 360;
    elseif (ang_temp(ii)-ang_temp(ii-1)) < -180
        ang_temp = ang_temp + 360;
    end

    insertion_angle.angle(ii) = ang_temp(ii);
end

insertion_angle.angle = -(insertion_angle.angle-360);

%% Smooth
insertion_angle.angle_smooth = smooth(insertion_angle.time, insertion_angle.angle, smooth_span,'rloess')';

%% Save
save(filename_save,'insertion_angle');

%% Plot
figure(); clf(); grid on; hold on;
xlabel('Time (s)')
ylabel('Insertion Angle (deg)')
plot(insertion_angle.time,insertion_angle.angle,'m','LineWidth',2)
plot(insertion_angle.time,insertion_angle.angle_smooth,'b','LineWidth',2)