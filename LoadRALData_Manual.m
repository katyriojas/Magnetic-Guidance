%% Import Data for Manual R-AL Insertion Trials
% Trevor Bruns and Katy Riojas
% Last Updated: 11/14/19

update_saved_manual_data = true; % specify whether to save regenerated data

%% Import Force data
addpath('functions','classDef','data');
load('avg_cal_slopes.mat'); % load force sensor calibration slopes
base_path = 'data'; % base path to csvs

filepaths_manual_phantom(1).force = fullfile(base_path, 'phantom\manual\phantom_manual_ea_trial1_force.csv');
filepaths_manual_phantom(2).force = fullfile(base_path, 'phantom\manual\phantom_manual_ea_trial2_force.csv');
filepaths_manual_phantom(3).force = fullfile(base_path, 'phantom\manual\phantom_manual_ea_trial3_force.csv');
filepaths_manual_phantom(4).force = fullfile(base_path, 'phantom\manual\phantom_manual_ea2_trial4_force.csv');

filepaths_manual_cadaver(1).force = fullfile(base_path, 'cadaver\manual\cadaver_manual_ea3_trial1_force.csv');
filepaths_manual_cadaver(2).force = fullfile(base_path, 'cadaver\manual\cadaver_manual_ea3_trial2_force.csv');
filepaths_manual_cadaver(3).force = fullfile(base_path, 'cadaver\manual\cadaver_manual_ea3_trial3_force.csv');


%% Create data objects
force_smooth_span = 40; % [# samples]

for ii = 1:length(filepaths_manual_phantom)
    data_manual_phantom(ii).nano = Nano17Data(filepaths_manual_phantom(ii).force, cal_slopes);
    data_manual_phantom(ii).nano.smooth_span = force_smooth_span;
end
for ii = 1:length(filepaths_manual_cadaver)
    data_manual_cadaver(ii).nano = Nano17Data(filepaths_manual_cadaver(ii).force, cal_slopes);
    data_manual_cadaver(ii).nano.smooth_span  = force_smooth_span;
end

%% Trim Trials
% Function "trimManualTrials" trims according to video footage
[phantom_trimmed,cadaver_trimmed,releaseTimes] =...
    trimManualTrials(data_manual_phantom,data_manual_cadaver);

%% Pull in phantom AID data (manually selected) 
% Note only have this for phantom
for ii = 1:length(filepaths_manual_phantom)
    load(strcat('data\phantom\manual\pman',num2str(ii),'_angular_depth.mat'));
    data_manual_phantom(ii).angular_depth = insertion_angle;
end

%% Interpolate force points at trimmed points
for ii = 1:length(filepaths_manual_phantom)
 
    data_manual_phantom(ii).interp_angdepth = ...
        interp1(data_manual_phantom(ii).angular_depth.time,...
                data_manual_phantom(ii).angular_depth.angle,...
                phantom_trimmed(ii).time - phantom_trimmed(ii).time(1),...
                'linear','extrap');
end

%% Add trimmed result to struct
for ii = 1:length(filepaths_manual_phantom)
    data_manual_phantom(ii).Fmag_trimmed = ...
        data_manual_phantom(ii).nano.Fmag(phantom_trimmed(ii).indices);
    data_manual_phantom(ii).time_trimmed = ...
        phantom_trimmed(ii).time - phantom_trimmed(ii).time(1);
    data_manual_phantom(ii).trim_idx = phantom_trimmed(ii).indices;
end

for ii = 1:length(filepaths_manual_cadaver)
    data_manual_cadaver(ii).Fmag_trimmed = ...
        data_manual_cadaver(ii).nano.Fmag(cadaver_trimmed(ii).indices);
    data_manual_cadaver(ii).time_trimmed = ...
        cadaver_trimmed(ii).time - cadaver_trimmed(ii).time(1);
    data_manual_cadaver(ii).trim_idx = cadaver_trimmed(ii).indices;
end

if update_saved_manual_data
    save('data\phantom\data_manual_phantom.mat','data_manual_phantom');
    save('data\cadaver\data_manual_cadaver.mat','data_manual_cadaver');
end