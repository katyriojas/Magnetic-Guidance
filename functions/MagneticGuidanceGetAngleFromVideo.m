%% MagneticGuidanceGetAngleFromVideo
%
%   insertion_angle.time  => [s] time since start of insertion
%   insertion_angle.angle => [deg] angle of the tip in cochlea frame
%   insertion_angle.angle_smooth => angle smoothed by smooth_span
%
% Trevor Bruns
% September 2019

function insertion_angle = MagneticGuidanceGetAngleFromVideo(video_path, smooth_span)

%% parameters

if nargin < 1
    % sample video location
    folderpath_vid = 'D:\Trevor\My Documents\MED lab\Cochlear R01\Mag Steering\Experiments\RAL\phantom_g_mea1_trial1_1.25';
    filename_vid = 'trial1-guided-mea1-1.25-tracked.MP4';
    video_path = fullfile(folderpath_vid, filename_vid);
end

if nargin < 2
    smooth_span = 0.02;
end

% reference marker colors
center.color = [1 0 0]; % red
angle0.color = [0 1 0]; % green
tip.color = [0 0 1]; % blue

% hue, saturation, value tolerances for segmenting colors
hsv_tolerance = [0.1, 0.3, 0.3];


%% find 0 degree angle reference

% load video and read first frame
vid = VideoReader(video_path);
vid.CurrentTime = 0;
curr_frame = readFrame(vid);

% center
center.pixels = segmentPixelsByColor(curr_frame, center.color, hsv_tolerance);
center.props = regionprops(center.pixels, 'basic');

% 0 degree marker
angle0.pixels = segmentPixelsByColor(curr_frame, angle0.color, hsv_tolerance);
angle0.props = regionprops(angle0.pixels, 'basic');

% ensure only 1 region was found
if (length(center.props)>1) || (length(angle0.props)>1)
    error('More than one matching region found, check your settings for ''hsv_tolerance''')
end

% 0 degree vector
angle0.vec = angle0.props.Centroid - center.props.Centroid;


%% step through video and determine angle
num_frames = floor(vid.FrameRate*vid.Duration) - 1; % approximate so we can pre-allocate without overshooting
insertion_angle.time  = zeros(1,num_frames); % [s]
insertion_angle.angle = zeros(1,num_frames); % [deg]
frame_count = 0;
vid.CurrentTime = 0;

% tic;
while hasFrame(vid)
    % load next frame
    curr_frame = readFrame(vid);
    frame_count = frame_count + 1;
    
    % segment tip marker
    tip.pixels = segmentPixelsByColor(curr_frame, tip.color, hsv_tolerance);
    tip.props = regionprops(tip.pixels, 'basic');
    
    % ensure only 1 region was found
    if (length(tip.props)>1)
        error('More than one matching region found, check your settings for ''hsv_tolerance''')
    end
    
    % find angle
    tip.vec = tip.props.Centroid - center.props.Centroid;
    insertion_angle.time(frame_count) = vid.CurrentTime;
    insertion_angle.angle(frame_count) = rad2deg( vectorAngle(angle0.vec, tip.vec) );
    
end
% toc

% ensure time increases monotonically (last frame shows previous time for some reason)
insertion_angle.time(end)  = [];
insertion_angle.angle(end) = [];

%% account 360 degree 'rollover'
ang_temp = insertion_angle.angle;
for ii=2:length(ang_temp)
    if (ang_temp(ii)-ang_temp(ii-1)) > 180
        ang_temp = ang_temp - 360;
    elseif (ang_temp(ii)-ang_temp(ii-1)) < -180
        ang_temp = ang_temp + 360;
    end

    insertion_angle.angle(ii) = ang_temp(ii);
end

%% smooth
insertion_angle.angle_smooth = smooth(insertion_angle.time, insertion_angle.angle, smooth_span, 'rloess')';

end

% %% plot
% figure(13); clf(13);
% plot(1.25*insertion_angle.time, insertion_angle.angle, 'm')
% hold on
% plot(1.25*insertion_angle.time, insertion_angle.angle_smooth, 'b')
% xlabel('time (s)')
% ylabel('insertion angle (deg)')