% Load in Cadaver Data
% Katy Riojas and Trevor Bruns
% Last Updated: 11/19/19

% This script will pull in the .csv files from the RA-L Robotic Cadaver 
% Experiments

% The outputs of the script will be data struct with each row as a 
% unguided (magnet-tipped array) /guided trio of trials 
% (script assumes equal number of trials for processing)

% A Magnetic_Guidance class instance is created for each trial

addpath('classDef','functions','data');
load('avg_cal_slopes.mat'); % loads cal_slopes
base_path = 'data\cadaver\';

filepaths_robotic_cadaver(1).nomag_mea.smaract = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial1_smaract.csv');
filepaths_robotic_cadaver(1).nomag_mea.force   = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial1_force.csv');

filepaths_robotic_cadaver(2).nomag_mea.smaract = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial2_smaract.csv');
filepaths_robotic_cadaver(2).nomag_mea.force   = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial2_force.csv');

filepaths_robotic_cadaver(3).nomag_mea.smaract = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial3_smaract.csv');
filepaths_robotic_cadaver(3).nomag_mea.force   = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial3_force.csv');

filepaths_robotic_cadaver(1).mag.smaract   = fullfile(base_path, 'guided\cadaver_g_mea2_trial1_smaract.csv');
filepaths_robotic_cadaver(1).mag.force     = fullfile(base_path, 'guided\cadaver_g_mea2_trial1_force.csv');

filepaths_robotic_cadaver(2).mag.smaract = fullfile(base_path, 'guided\cadaver_g_mea2_trial2_smaract.csv');
filepaths_robotic_cadaver(2).mag.force   = fullfile(base_path, 'guided\cadaver_g_mea2_trial2_force.csv');

filepaths_robotic_cadaver(3).mag.smaract = fullfile(base_path, 'guided\cadaver_g_mea2_trial3_smaract.csv');
filepaths_robotic_cadaver(3).mag.force   = fullfile(base_path, 'guided\cadaver_g_mea2_trial3_force.csv');

%% Create data objects
force_smooth_span = 40; % set smooth span [# samples]

for ii = 1:length(filepaths_robotic_cadaver)
    data_robotic_cadaver(ii).nomag = MagneticGuidanceData(filepaths_robotic_cadaver(ii).nomag_mea, cal_slopes);
    data_robotic_cadaver(ii).mag = MagneticGuidanceData(filepaths_robotic_cadaver(ii).mag, cal_slopes);
    data_robotic_cadaver(ii).nomag.smooth_span = force_smooth_span;
    data_robotic_cadaver(ii).mag.smooth_span   = force_smooth_span;
end

%% Trim to 125mN Max Force because this was the max allowable force, but
% there are still some data points post this force
robotic_cadaver_nomag_end = zeros(1,length(filepaths_robotic_cadaver)); % initialize a vector with end indices for nomag trials
robotic_cadaver_mag_end = zeros(1,length(filepaths_robotic_cadaver)); % initialize a vector with end indices for mag trials

Fthresh = 125; %[mN] force threshold

% Just pulling out end index - note, this won't actually trim the data
for ii = 1:length(filepaths_robotic_cadaver)
    nomagii_end = find(data_robotic_cadaver(ii).nomag.Fmag>Fthresh,1);
    magii_end = find(data_robotic_cadaver(ii).mag.Fmag>Fthresh,1);
    if isempty(nomagii_end)
        robotic_cadaver_nomag_end(ii) = length(data_robotic_cadaver(ii).nomag.Fmag);
    else
        robotic_cadaver_nomag_end(ii) = nomagii_end;
    end
    if isempty(magii_end)
        robotic_cadaver_mag_end(ii) = length(data_robotic_cadaver(ii).mag.Fmag);
    else
        robotic_cadaver_mag_end(ii) = magii_end;
    end
end