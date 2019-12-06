% Process RA-L Phantom Data

% This script will pull in the .csv files from the RA-L Robotic
% Phantom Experiments and postprocessed .mp4 files

% The outputs of the script will be data struct with each row featuring 
% a different trial (4 total) and the columns specifying method of
% insertion, that is, unguided (non-magnet tipped array)/unguided (MEA)/guided
% (script assumes equal number of trials for processing)

% A Magnetic_Guidance class instance is created for each trial
% The videos are processed to analyze AID vs. linear insertion depth

update_saved_phantom_data = true; % default is true - update saved.mat after regenerating


%%
addpath('classDef','functions','data');
load('avg_cal_slopes.mat'); % loads cal_slopes
base_path = 'data\phantom\';

% Unguided, unmodified EA
filepaths_robotic_phantom(1).nomag_ea.smaract = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial1_1.25_smaract.csv');
filepaths_robotic_phantom(1).nomag_ea.force   = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial1_1.25_force.csv');

filepaths_robotic_phantom(2).nomag_ea.smaract = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial2_1.25_smaract.csv');
filepaths_robotic_phantom(2).nomag_ea.force   = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial2_1.25_force.csv');

filepaths_robotic_phantom(3).nomag_ea.smaract = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial3_1.25_smaract.csv');
filepaths_robotic_phantom(3).nomag_ea.force   = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial3_1.25_force.csv');

filepaths_robotic_phantom(4).nomag_ea.smaract = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial4_1.25_smaract.csv');
filepaths_robotic_phantom(4).nomag_ea.force   = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial4_1.25_force.csv');

% Unguided, magnet-tipped EA
filepaths_robotic_phantom(1).nomag_mea.smaract = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial1_1.25_smaract.csv');
filepaths_robotic_phantom(1).nomag_mea.force   = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial1_1.25_force.csv');

filepaths_robotic_phantom(2).nomag_mea.smaract = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial2_1.25_smaract.csv');
filepaths_robotic_phantom(2).nomag_mea.force   = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial2_1.25_force.csv');

filepaths_robotic_phantom(3).nomag_mea.smaract = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial3_1.25_smaract.csv');
filepaths_robotic_phantom(3).nomag_mea.force   = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial3_1.25_force.csv');

filepaths_robotic_phantom(4).nomag_mea.smaract = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial4_1.25_smaract.csv');
filepaths_robotic_phantom(4).nomag_mea.force   = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial4_1.25_force.csv');

% Guided (with magnet)
filepaths_robotic_phantom(1).mag.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial1_1.25_smaract.csv');
filepaths_robotic_phantom(1).mag.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial1_1.25_force.csv');

filepaths_robotic_phantom(2).mag.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial2_1.25_smaract.csv');
filepaths_robotic_phantom(2).mag.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial2_1.25_force.csv');

filepaths_robotic_phantom(3).mag.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial3_1.25_smaract.csv');
filepaths_robotic_phantom(3).mag.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial3_1.25_force.csv');

filepaths_robotic_phantom(4).mag.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial4_1.25_smaract.csv');
filepaths_robotic_phantom(4).mag.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial4_1.25_force.csv');


%% Create data objects
force_smooth_span = 40; % set smooth span[# samples]

for ii=1:length(filepaths_robotic_phantom)
    data_robotic_phantom(ii).nomag_ea  = MagneticGuidanceData(filepaths_robotic_phantom(ii).nomag_ea,  cal_slopes);
    data_robotic_phantom(ii).nomag_mea = MagneticGuidanceData(filepaths_robotic_phantom(ii).nomag_mea, cal_slopes);
    data_robotic_phantom(ii).mag       = MagneticGuidanceData(filepaths_robotic_phantom(ii).mag,       cal_slopes);
    data_robotic_phantom(ii).nomag_ea.smooth_span  = force_smooth_span;
    data_robotic_phantom(ii).nomag_mea.smooth_span = force_smooth_span;
    data_robotic_phantom(ii).mag.smooth_span       = force_smooth_span;
end


%% Determine angular insertion depths from processed video

angle_smooth_span = 20;

% video_base_path = 'C:\Users\riojaske\Documents\magsteer\Videos';
video_base_path = 'D:\Trevor\My Documents\MED lab\Cochlear R01\Mag Steering\Experiments\RAL';

filepaths_robotic_phantom(1).nomag_mea.vid = fullfile(video_base_path, 'phantom_ug_mea1_trial1_1.25\trial1-ug-mea1-tracked.mp4');
filepaths_robotic_phantom(1).mag.vid       = fullfile(video_base_path, 'phantom_g_mea1_trial1_1.25\trial1-guided-mea1-1.25-tracked.mp4');

filepaths_robotic_phantom(2).nomag_mea.vid = fullfile(video_base_path, 'phantom_ug_mea1_trial2_1.25\trial2-ug-mea1-tracked.mp4');
filepaths_robotic_phantom(2).mag.vid       = fullfile(video_base_path, 'phantom_g_mea1_trial2_1.25\trial2-g-mea1-tracked.mp4');

filepaths_robotic_phantom(3).nomag_mea.vid = fullfile(video_base_path, 'phantom_ug_mea1_trial3_1.25\trial3-ug-mea1-tracked.mp4');
filepaths_robotic_phantom(3).mag.vid       = fullfile(video_base_path, 'phantom_g_mea1_trial3_1.25\trial3-g-mea1-tracked.mp4');

filepaths_robotic_phantom(4).nomag_mea.vid = fullfile(video_base_path, 'phantom_ug_mea1_trial4_1.25\trial4-ug-mea1-tracked.mp4');
filepaths_robotic_phantom(4).mag.vid       = fullfile(video_base_path, 'phantom_g_mea1_trial4_1.25\trial4-g-mea1-tracked.mp4');

for ii = 1:length(filepaths_robotic_phantom)
    
    % Segment video frames to determine angular depth
    data_robotic_phantom(ii).mag_angular_depth = MagneticGuidanceGetAngleFromVideo(filepaths_robotic_phantom(ii).mag.vid, angle_smooth_span);
    data_robotic_phantom(ii).nomag_mea_angular_depth = MagneticGuidanceGetAngleFromVideo(filepaths_robotic_phantom(ii).nomag_mea.vid, angle_smooth_span);
    
    % Interpolate to find angular depth at each force measurement time
    data_robotic_phantom(ii).mag_interp_angdepth = ...
        interp1(data_robotic_phantom(ii).mag_angular_depth.time,...
                data_robotic_phantom(ii).mag_angular_depth.angle_smooth,...
                data_robotic_phantom(ii).mag.time_insertion);

    if ii~=3 % video stopped early on nomag trial 3  
        data_robotic_phantom(ii).nomag_mea_interp_angdepth = ...
            interp1(data_robotic_phantom(ii).nomag_mea_angular_depth.time,...
                    data_robotic_phantom(ii).nomag_mea_angular_depth.angle_smooth,...
                    data_robotic_phantom(ii).nomag_mea.time_insertion);
    end

end


%% Use nomag trial 4 to fill remainder of nomag trial 3

% trial 3 data until stop point
nomag3_vidstop_time  = data_robotic_phantom(3).nomag_mea_angular_depth.time(end);
nomag3_stop_ind = find(data_robotic_phantom(3).nomag_mea.time_insertion > nomag3_vidstop_time, 1);

data_robotic_phantom(3).nomag_mea_interp_angdepth = interp1(data_robotic_phantom(3).nomag_mea_angular_depth.time,...
                                            data_robotic_phantom(3).nomag_mea_angular_depth.angle_smooth,...
                                            data_robotic_phantom(3).nomag_mea.time_insertion(1:nomag3_stop_ind), 'linear', 'extrap');

% addition from trial 4
nomag3_vidstop_depth = data_robotic_phantom(3).nomag_mea.depth_insertion( find(data_robotic_phantom(3).nomag_mea.time_insertion > nomag3_vidstop_time, 1) );
nomag3_total_time = data_robotic_phantom(3).nomag_mea.time_insertion(end);
nomag3_total_depth = data_robotic_phantom(3).nomag_mea.depth_insertion(end);


nomag4_inds = find(data_robotic_phantom(4).nomag_mea.depth_insertion > nomag3_vidstop_depth, 1): find(data_robotic_phantom(4).nomag_mea.depth_insertion > nomag3_total_depth, 1);

nomag3_interp_angdepth = [data_robotic_phantom(3).nomag_mea_interp_angdepth; data_robotic_phantom(4).nomag_mea_interp_angdepth(nomag4_inds)];

nomag3_insertion_depths = linspace(data_robotic_phantom(3).nomag_mea.depth_insertion(1), data_robotic_phantom(3).nomag_mea.depth_insertion(end), length(nomag3_interp_angdepth));

% Re-interpolate so you have angular insertion depths match linear
% Insertion depths from ROS
data_robotic_phantom(3).nomag_mea_interp_angdepth = ...
    interp1(nomag3_insertion_depths,nomag3_interp_angdepth,...
            data_robotic_phantom(3).nomag_mea.depth_insertion, 'linear', 'extrap');


%% Binning

degspan = 2; % [deg] span/width of each bin

% compute the largest angular insertion depth to bound our bins
max_ang = 0;
for ii = 1:length(data_robotic_phantom)
    max_ang = max([max_ang,...
                   max(data_robotic_phantom(ii).nomag_mea_interp_angdepth),...
                   max(data_robotic_phantom(ii).mag_interp_angdepth)]);
end
num_bins = ceil(max_ang/degspan);           

% create vector of the bin edges
bin_edges = 0:degspan:num_bins*degspan;

% create vector of the bin centers
bins = bin_edges(2:end) - degspan/2;

% sort into bins
for ii = 1:length(data_robotic_phantom)

    % determine the bin for each interp_angdepth point
    [~,~,data_robotic_phantom(ii).nomag_mea_binned.ind] = histcounts(data_robotic_phantom(ii).nomag_mea_interp_angdepth, bin_edges);
    [~,~,data_robotic_phantom(ii).mag_binned.ind]       = histcounts(data_robotic_phantom(ii).mag_interp_angdepth, bin_edges);

    % compute mean of all the force measurements within each bin (NaN for bins with no measurements)
    data_robotic_phantom(ii).nomag_mea_binned.Fmean = zeros(size(bins));
    data_robotic_phantom(ii).mag_binned.Fmean = zeros(size(bins));
    for jj = 1:length(bins)
        data_robotic_phantom(ii).nomag_mea_binned.Fmean(jj) = ...
            mean( data_robotic_phantom(ii).nomag_mea.Fmag( data_robotic_phantom(ii).nomag_mea_binned.ind == jj ) );

        data_robotic_phantom(ii).mag_binned.Fmean(jj) = ...
            mean( data_robotic_phantom(ii).mag.Fmag( data_robotic_phantom(ii).mag_binned.ind == jj ) );
    end
    
    % also save the actual bins and edges
    data_robotic_phantom(ii).nomag_mea_binned.bins  = bins;
    data_robotic_phantom(ii).nomag_mea_binned.edges = bin_edges;
    data_robotic_phantom(ii).mag_binned.bins  = bins;
    data_robotic_phantom(ii).mag_binned.edges = bin_edges;
end



%% Update Saved Phantom Data if it is called for
if update_saved_phantom_data
   save('data\phantom\data_robotic_phantom.mat','data_robotic_phantom');
end