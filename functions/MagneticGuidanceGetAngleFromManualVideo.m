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

% filename_vid = '\videos\trial1\trial1-manual-ea1-cropped-10fps.mp4';
% filename_save = strcat(folderpath_vid,'\pman1_angular_depth.mat');
% filename_vid = '\videos\trial2\trial2-manual-ea1-10fps.mp4';
% filename_save = strcat(folderpath_vid,'\pman2_angular_depth.mat');
% filename_vid = '\videos\trial3\trial3-manual-ea1-cropped-10fps.mp4';
% filename_save = strcat(folderpath_vid,'\pman3_angular_depth.mat');
filename_vid = '\videos\trial4\trial4-manual-ea2-cropped-10fps.mp4';
filename_save = strcat(folderpath_vid,'\pman4_angular_depth.mat');

video_path = fullfile(folderpath_vid, filename_vid);
smooth_span = 10;
magnification = 300;

% load video and read first frame
vid = VideoReader(video_path);
vid.CurrentTime = 0;
curr_frame = readFrame(vid);
imshow(curr_frame,'InitialMagnification',magnification); % show frame

% center
fprintf('Select Center Location\n');
[centerX,centerY] = ginput(1);
centerPnt = [centerX;centerY];

fprintf('Select 0 degree Mark\n');
% 0 degree marker
[angle0X,angle0Y] = ginput(1);
angle0Pnt = [angle0X;angle0Y];
clf();

% 0 degree vector
angle0.vec = angle0Pnt - centerPnt;

%% Step through video and determine angle
num_frames = floor(vid.FrameRate*vid.Duration) - 1; % approximate so we can pre-allocate without overshooting
insertion_angle.time  = zeros(1,num_frames); % [s]
insertion_angle.angle = zeros(1,num_frames); % [deg]
frame_count = 0;
vid.CurrentTime = 0;

% tic;
while hasFrame(vid)
    if frame_count == 0
      fprintf('Select Tip location (most distal end of most apical electrode pad\n');
    end
    % load next frame
    curr_frame = readFrame(vid);
    imshow(curr_frame,'InitialMagnification',magnification); % show frame
    
    frame_count = frame_count + 1;
    percDone = 100*frame_count/num_frames;
    fprintf('%2.1f perc\n',percDone)
    
    % segment tip marker
    [tipX,tipY] = ginput(1);
    tipPnt = [tipX;tipY];
    
    % find angle
    tip.vec = tipPnt - centerPnt;
    insertion_angle.time(frame_count) = vid.CurrentTime;
    insertion_angle.angle(frame_count) = rad2deg(vectorAngle(angle0.vec',tip.vec'));
    
end
% toc

% ensure time increases monotonically (last frame shows previous time for some reason)
insertion_angle.time(end)  = [];
insertion_angle.angle(end) = [];

%% account 360 degree 'rollover'
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

%% smooth
insertion_angle.angle_smooth = smooth(insertion_angle.time, insertion_angle.angle, smooth_span, 'rloess')';

%% save
save(filename_save,'insertion_angle');

%% plot
figure(); clf(); grid on; hold on;
xlabel('time (s)')
ylabel('insertion angle (deg)')
plot(insertion_angle.time, insertion_angle.angle, 'm')
plot(insertion_angle.time, insertion_angle.angle_smooth, 'b')