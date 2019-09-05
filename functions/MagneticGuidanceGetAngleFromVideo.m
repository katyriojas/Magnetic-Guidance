%% MagneticGuidanceGetAngleFromVideo
%
% Trevor Bruns
% September 2019

%% user parameters

% video location
folderpath_vid = 'D:\Trevor\My Documents\MED lab\Cochlear R01\Mag Steering\Experiments\RAL\phantom_g_mea1_trial1_1.25';
filename_vid     = 'trial1-guided-mea1-1.25-tracked.MP4';

% CIE L*a*b threshold values for center (red), 0 degree mark (green), tip (blue)
center.a = 55;
center.b = 25;
angle0.a = 

%% load video

vid = VideoReader(fullfile(folderpath_vid, filename_vid));
frame_step = 1/vid.FrameRate; % [s] time step between frames
vid.CurrentTime = 0;

% figure(12); clf(12);
% imshow(curr_frame, 'Border','tight');
% hold on


%%


curr_frame = readFrame(vid);
