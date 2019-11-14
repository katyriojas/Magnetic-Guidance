% Katy Riojas
% Plot angle and video
% Adjust to save center, zero, and xy positions

% Sample video location
clear all; close all; clc;
folderpath_vid = 'C:\Users\riojaske\Documents\magsteer\Magnetic-Guidance\data\phantom\manual';

filename_vid = '\videos\trial1\trial1-manual-ea1-cropped-10fps.mp4';
filename_save = strcat(folderpath_vid,'\pman1_angular_depth.mat');
% filename_vid = '\videos\trial2\trial2-manual-ea1-10fps.mp4';
% filename_save = strcat(folderpath_vid,'\pman2_angular_depth.mat');
% filename_vid = '\videos\trial3\trial3-manual-ea1-cropped-10fps.mp4';
% filename_save = strcat(folderpath_vid,'\pman3_angular_depth.mat');
% filename_vid = '\videos\trial4\trial4-manual-ea2-cropped-10fps.mp4';
% filename_save = strcat(folderpath_vid,'\pman4_angular_depth.mat');
load(filename_save);

mag = 300;

video_path = fullfile(folderpath_vid, filename_vid);
% Load video and read first frame
vid = VideoReader(video_path);
vid.CurrentTime = 0;
curr_frame = readFrame(vid);
imshow(curr_frame,'InitialMagnification',mag); % show frame

% Center
fprintf('Select Center Location\n');
[centerX,centerY] = ginput(1);
centerPnt = [centerX;centerY];
fprintf('Select 0 degree Mark\n');
% 0 degree marker
[angle0X,angle0Y] = ginput(1);
angle0Pnt = [angle0X;angle0Y];
clf();

% 0 degree vector
angle0vec = angle0Pnt - centerPnt;

num_frames = floor(vid.FrameRate*vid.Duration) - 1; % approximate so we can pre-allocate without overshooting
% insertion_angle.time  = zeros(1,num_frames); % [s]
% insertion_angle.angle = zeros(1,num_frames); % [deg]
frame_count = 0;
vid.CurrentTime = 0;

% tic;
while hasFrame(vid)

    % load next frame
    curr_frame = readFrame(vid);
    imshow(curr_frame,'InitialMagnification',mag); % show frame
    hold on;
    frame_count = frame_count + 1;
    percDone = 100*frame_count/num_frames;
    fprintf('%2.1f perc\n',percDone)
    
    % Plot angle - first get vector
    line_length = 100;
    Rz = roz(deg2rad(-insertion_angle.angle(frame_count)));
    zeroVec = angle0vec + centerPnt;
    h = Rz(1:2,1:2)*angle0vec + centerPnt;
    
    %Plot
    line([centerPnt(1),h(1)],[centerPnt(2),h(2)]);
    line([centerPnt(1),zeroVec(1)],[centerPnt(2),zeroVec(2)]);
    scatter(centerPnt(1),centerPnt(2),'filled','b');
    scatter(h(1),h(2),'filled','r');
    scatter(zeroVec(1),zeroVec(2),'filled','g');
    text(381,144,num2str(insertion_angle.angle(frame_count)));
    text(381,160,num2str(insertion_angle.angle_smooth(frame_count)));
    drawnow;
    clf();
    
end

figure(); grid on; hold on;
plot(insertion_angle.angle,'b:','LineWidth',1)
plot(insertion_angle.angle_smooth,'r:','LineWidth',1)
legend('Angle','Smoothed Angle');
legend('Raw Angle','Smoothed Angle');