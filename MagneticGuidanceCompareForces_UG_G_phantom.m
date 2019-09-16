%% Data Analysis for Phantom Trials
% Trevor Bruns and Katy Riojas
% Last Updated: 9/5/19
clear all; clc;

%% INPUTS
addpath('classDef','functions','data');
load('avg_cal_slopes.mat');
load('data\phantom\phantom_T_st_fixture.mat');
T_fixture_st = inv(T_st_fixture);
smooth_span = 0.06; % specify smooth span

% Filepaths to CSVs exported from ROS bags
base_path = 'data\phantom\';

filepaths_nomag1_mea.smaract = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial1_1.25_smaract.csv');
filepaths_nomag1_mea.force   = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial1_1.25_force.csv');

filepaths_nomag2_mea.smaract = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial2_1.25_smaract.csv');
filepaths_nomag2_mea.force   = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial2_1.25_force.csv');

filepaths_nomag3_mea.smaract = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial3_1.25_smaract.csv');
filepaths_nomag3_mea.force   = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial3_1.25_force.csv');

filepaths_nomag4_mea.smaract = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial4_1.25_smaract.csv');
filepaths_nomag4_mea.force   = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial4_1.25_force.csv');

filepaths_mag1.smaract   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial1_1.25_smaract.csv');
filepaths_mag1.force     = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial1_1.25_force.csv');

filepaths_mag2.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial2_1.25_smaract.csv');
filepaths_mag2.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial2_1.25_force.csv');

filepaths_mag3.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial3_1.25_smaract.csv');
filepaths_mag3.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial3_1.25_force.csv');

filepaths_mag4.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial4_1.25_smaract.csv');
filepaths_mag4.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial4_1.25_force.csv');

filepaths_nomag1_ea.smaract = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial1_1.25_smaract.csv');
filepaths_nomag1_ea.force   = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial1_1.25_force.csv');

filepaths_nomag2_ea.smaract = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial2_1.25_smaract.csv');
filepaths_nomag2_ea.force   = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial2_1.25_force.csv');

filepaths_nomag3_ea.smaract = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial3_1.25_smaract.csv');
filepaths_nomag3_ea.force   = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial3_1.25_force.csv');

filepaths_nomag4_ea.smaract = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial4_1.25_smaract.csv');
filepaths_nomag4_ea.force   = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial4_1.25_force.csv');


%% create data objects
data_nomag1_mea = MagneticGuidanceData(filepaths_nomag1_mea, cal_slopes);
data_nomag2_mea = MagneticGuidanceData(filepaths_nomag2_mea, cal_slopes);
data_nomag3_mea = MagneticGuidanceData(filepaths_nomag3_mea, cal_slopes);
data_nomag4_mea = MagneticGuidanceData(filepaths_nomag4_mea, cal_slopes);

data_mag1 = MagneticGuidanceData(filepaths_mag1, cal_slopes);
data_mag2 = MagneticGuidanceData(filepaths_mag2, cal_slopes);
data_mag3 = MagneticGuidanceData(filepaths_mag3, cal_slopes);
data_mag4 = MagneticGuidanceData(filepaths_mag4, cal_slopes);

data_nomag1_ea = MagneticGuidanceData(filepaths_nomag1_ea, cal_slopes);
data_nomag2_ea = MagneticGuidanceData(filepaths_nomag2_ea, cal_slopes);
data_nomag3_ea = MagneticGuidanceData(filepaths_nomag3_ea, cal_slopes);
data_nomag4_ea = MagneticGuidanceData(filepaths_nomag4_ea, cal_slopes);


%% set smoothing span (proportion of data points)
force_smooth_span = 40; % [# samples]

data_nomag1_mea.smooth_span = force_smooth_span;
data_nomag2_mea.smooth_span = force_smooth_span;
data_nomag3_mea.smooth_span = force_smooth_span;
data_nomag4_mea.smooth_span = force_smooth_span;

data_mag1.smooth_span = force_smooth_span;
data_mag2.smooth_span = force_smooth_span;
data_mag3.smooth_span = force_smooth_span;
data_mag4.smooth_span = force_smooth_span;

data_nomag1_ea.smooth_span = force_smooth_span;
data_nomag2_ea.smooth_span = force_smooth_span;
data_nomag3_ea.smooth_span = force_smooth_span;
data_nomag4_ea.smooth_span = force_smooth_span;


%% determine angular insertion depths from processed video
% angle_smooth_span = 20;
% mag1.angular_depth = MagneticGuidanceGetAngleFromVideo(...
%     'D:\Trevor\My Documents\MED lab\Cochlear R01\Mag Steering\Experiments\RAL\phantom_g_mea1_trial1_1.25\trial1-guided-mea1-1.25-tracked.MP4'...
%     ,angle_smooth_span);
% 
% nomag1.angular_depth = MagneticGuidanceGetAngleFromVideo(...
%     'D:\Trevor\My Documents\MED lab\Cochlear R01\Mag Steering\Experiments\RAL\phantom_ug_mea1_trial1_1.25\trial1-ug-mea1-tracked.mp4'...
%     ,angle_smooth_span);
% 
% % interpolate to find angular depth at each force measurement time
% mag1.interp_angdepth   = interp1(mag1.angular_depth.time, mag1.angular_depth.angle, data_mag1.time_insertion, 'linear', 'extrap');
% nomag1.interp_angdepth = interp1(nomag1.angular_depth.time, nomag1.angular_depth.angle, data_nomag1_mea.time_insertion, 'linear', 'extrap');


%% testing
% force_smooth_span = 40;
% data_mag1.smooth_span = force_smooth_span;
% 
% figure(1); clf(1); 
% hold on; grid on;
% 
% plot(mag1.interp_angdepth, data_mag1.Fx,          'Color', [1,0,0, 0.3], 'LineWidth', 1);
% h_mag1(1) = plot(mag1.interp_angdepth, data_mag1.Fx_smooth,   'Color', [1,0,0, 1.0], 'LineWidth', 2);
% plot(mag1.interp_angdepth, data_mag1.Fy,          'Color', [0,1,0, 0.3], 'LineWidth', 1);
% h_mag1(2) = plot(mag1.interp_angdepth, data_mag1.Fy_smooth,   'Color', [0,1,0, 1.0], 'LineWidth', 2);
% plot(mag1.interp_angdepth, data_mag1.Fz,          'Color', [0,0,1, 0.3], 'LineWidth', 1);
% h_mag1(3) = plot(mag1.interp_angdepth, data_mag1.Fz_smooth,   'Color', [0,0,1, 1.0], 'LineWidth', 2);

% figure; hold on
% plot(nomag1.angular_depth.time, nomag1.angular_depth.angle, 'r')
% plot(nomag1.angular_depth.time, nomag1.angular_depth.angle_smooth, 'b')

%%
% figure(1); clf(1); 
% hold on; grid on;
% 
% plot(mag1.interp_angdepth, data_mag1.Fx,          'Color', [1,0,0, 0.3], 'LineWidth', 1);
% h_mag1(1) = plot(mag1.interp_angdepth, data_mag1.Fx_smooth,   'Color', [1,0,0, 1.0], 'LineWidth', 2);
% plot(mag1.interp_angdepth, data_mag1.Fy,          'Color', [0,1,0, 0.3], 'LineWidth', 1);
% h_mag1(2) = plot(mag1.interp_angdepth, data_mag1.Fy_smooth,   'Color', [0,1,0, 1.0], 'LineWidth', 2);
% plot(mag1.interp_angdepth, data_mag1.Fz,          'Color', [0,0,1, 0.3], 'LineWidth', 1);
% h_mag1(3) = plot(mag1.interp_angdepth, data_mag1.Fz_smooth,   'Color', [0,0,1, 1.0], 'LineWidth', 2);
% 
% plot(nomag1.interp_angdepth, data_nomag1_mea.Fx,          'Color', [1,0,0, 0.3], 'LineWidth', 1);
% h_mag1(1) = plot(nomag1.interp_angdepth, data_nomag1_mea.Fx_smooth,   'Color', [1,0,0, 0.3], 'LineWidth', 2);
% plot(nomag1.interp_angdepth, data_nomag1_mea.Fy,          'Color', [0,1,0, 0.3], 'LineWidth', 1);
% h_mag1(2) = plot(nomag1.interp_angdepth, data_nomag1_mea.Fy_smooth,   'Color', [0,1,0, 0.3], 'LineWidth', 2);
% plot(nomag1.interp_angdepth, data_nomag1_mea.Fz,          'Color', [0,0,1, 0.3], 'LineWidth', 1);
% h_mag1(3) = plot(nomag1.interp_angdepth, data_nomag1_mea.Fz_smooth,   'Color', [0,0,1, 0.3], 'LineWidth', 2);
% 
% title('Force vs Angular Insertion Depth')
% xlabel('angular insertion depth (deg)')
% ylabel('force (mN)')
% legend(h_mag1, {'Fx', 'Fy', 'Fz'}, 'Location', 'sw')


% figure(2); clf(2);
% hold on; grid on;
% plot(mag1.interp_angdepth, data_mag1.Fmag,        'Color', [0,0,1, 0.3], 'LineWidth', 1);
% plot(mag1.interp_angdepth, data_mag1.Fmag_smooth, 'Color', [0,0,1, 1.0], 'LineWidth', 2);
% plot(nomag1.interp_angdepth, data_nomag1_mea.Fmag,        'Color', [0,0,0, 0.3], 'LineWidth', 1);
% plot(nomag1.interp_angdepth, data_nomag1_mea.Fmag_smooth, 'Color', [0,0,0, 1.0], 'LineWidth', 2);
% title('Force vs Angular Insertion Depth')
% xlabel('angular insertion depth (deg)')
% ylabel('force (mN)')


% subplot(2,1,2);
% hold on; grid on;
% % plot(data_mag1.depth_insertion, data_mag1.Fmag,        'Color', [0,0,0, 0.3], 'LineWidth', 1);
% % plot(data_mag1.depth_insertion, data_mag1.Fmag_smooth, 'Color', [0,0,0, 1.0], 'LineWidth', 2);
% plot(data_mag1.depth_insertion, data_mag1.Fx,          'Color', [1,0,0, 0.3], 'LineWidth', 1);
% h_mag1(1) = plot(data_mag1.depth_insertion, data_mag1.Fx_smooth,   'Color', [1,0,0, 1.0], 'LineWidth', 2);
% plot(data_mag1.depth_insertion, data_mag1.Fy,          'Color', [0,1,0, 0.3], 'LineWidth', 1);
% h_mag1(2) = plot(data_mag1.depth_insertion, data_mag1.Fy_smooth,   'Color', [0,1,0, 1.0], 'LineWidth', 2);
% plot(data_mag1.depth_insertion, data_mag1.Fz,          'Color', [0,0,1, 0.3], 'LineWidth', 1);
% h_mag1(3) = plot(data_mag1.depth_insertion, data_mag1.Fz_smooth,   'Color', [0,0,1, 1.0], 'LineWidth', 2);
% title('Force vs Linear Insertion Depth')
% xlabel('linear insertion depth (mm)')
% ylabel('force (mN)')
% legend(h_mag1, {'Fx', 'Fy', 'Fz'}, 'Location', 'sw')

%% Trim data to when the array is within 10 degrees of final placement
nomag1_tstart_vid = 2.12;
nomag1_tend_vid = 20.53;
nomag1_tend = nomag1_tend_vid - nomag1_tstart_vid;

nomag2_tstart_vid = 4.48;
nomag2_tend_vid = 24.43;
nomag2_tend = nomag2_tend_vid - nomag2_tstart_vid;

% Video cut off bc this was the trial where the magnet shut off, 
%just taking the end time as the end of insertion for this trial
nomag3_tstart_vid = 4.28;
nomag3_tend_vid = data_nomag3_mea.time_insertion(end);
nomag3_tend = nomag3_tend_vid - nomag3_tstart_vid;

nomag4_tstart_vid = 3.22;
nomag4_tend_vid = 24.25;
nomag4_tend = nomag4_tend_vid - nomag4_tstart_vid;

mag1_tstart_vid = 2.87;
mag1_tend_vid = 21.77;
mag1_tend = mag1_tend_vid - mag1_tstart_vid;

mag2_tstart_vid = 6.03;
mag2_tend_vid = 25.15;
mag2_tend = mag2_tend_vid - mag2_tstart_vid;

mag3_tstart_vid = 6.17;
mag3_tend_vid = 27.42;
mag3_tend =  mag3_tend_vid - mag3_tstart_vid;

mag4_tstart_vid = 4.05;
mag4_tend_vid = 25.3;
mag4_tend = mag4_tend_vid - mag4_tstart_vid;

nomag1_endi = find(data_nomag1_mea.time_insertion > nomag1_tend,1);
if isempty(nomag1_endi), nomag1_endi = length(data_nomag1_mea.time_insertion); end
nomag2_endi = find(data_nomag2_mea.time_insertion > nomag2_tend,1);
if isempty(nomag2_endi), nomag2_endi = length(data_nomag2_mea.time_insertion); end
nomag3_endi = find(data_nomag3_mea.time_insertion > nomag3_tend,1);
if isempty(nomag3_endi), nomag3_endi = length(data_nomag3_mea.time_insertion); end
nomag4_endi = find(data_nomag4_mea.time_insertion > nomag4_tend,1);
if isempty(nomag4_endi), nomag4_endi = length(data_nomag4_mea.time_insertion); end

mag1_endi = find(data_mag1.time_insertion > mag1_tend,1);
if isempty(mag1_endi), mag1_endi = length(data_mag1.time_insertion); end
mag2_endi = find(data_mag2.time_insertion > mag2_tend,1);
if isempty(mag2_endi), mag2_endi = length(data_mag2.time_insertion); end
mag3_endi = find(data_mag3.time_insertion > mag3_tend,1);
if isempty(mag3_endi), mag3_endi = length(data_mag3.time_insertion); end
mag4_endi = find(data_mag4.time_insertion > mag4_tend,1);
if isempty(mag4_endi), mag4_endi = length(data_mag4.time_insertion); end

%% Plot Props
alpha = 1; % line opacity
colorsMat = distinguishable_colors(12);
xyzColor = distinguishable_colors(3); % xyz colors
plotXYZData = 0; %toggle whether to plot XYZ data
maxY = max([data_nomag1_mea.Fmag;data_nomag2_mea.Fmag;...
            data_nomag3_mea.Fmag;data_nomag4_mea.Fmag;data_mag1.Fmag;...
            data_mag2.Fmag;data_mag3.Fmag;data_mag4.Fmag]);

%% Figure 1: Force Magnitude Trial Comparison
figure(1); clf(1); hold on; grid on;
title('Guided vs Unguided Insertion');
xlabel('Insertion Depth [mm]');
ylabel('||Force|| [mN]');
figure(1); clf(1);
subplot(2,4,1); grid on; hold on;
sgtitle('Raw and Smoothed Fmag Phantom Trials');
ylabel('Unguided ||Force|| [mN]');
ylim([0 maxY]);
% xlabel('Insertion Depth [mm]'); 
title('Trial 1');
plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fmag, 'Color', [colorsMat(1,:),  0.3*alpha], 'LineWidth',1);
plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth',1);
scatter(data_nomag1_mea.depth_insertion(nomag1_endi),data_nomag1_mea.Fmag(nomag1_endi),10,'filled','r');

subplot(2,4,2); grid on; hold on;
% xlabel('Insertion Depth [mm]'); 
title('Trial 2');
ylim([0 maxY]);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fmag, 'Color', [colorsMat(1,:), 0.3*alpha], 'LineWidth',1);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth',1);
scatter(data_nomag2_mea.depth_insertion(nomag2_endi),data_nomag2_mea.Fmag(nomag2_endi),10,'filled','r');

subplot(2,4,3); grid on; hold on;
% xlabel('Insertion Depth [mm]'); 
title('Trial 3');
ylim([0 maxY]);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fmag, 'Color', [colorsMat(1,:),  0.3*alpha], 'LineWidth',1);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth',1);
scatter(data_nomag3_mea.depth_insertion(nomag3_endi),data_nomag3_mea.Fmag(nomag3_endi),10,'filled','r');

subplot(2,4,4); grid on; hold on;
% xlabel('Insertion Depth [mm]'); 
title('Trial 4');
ylim([0 maxY]);
plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fmag, 'Color', [colorsMat(1,:),  0.3*alpha], 'LineWidth',1);
plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth',1);
scatter(data_nomag4_mea.depth_insertion(nomag4_endi),data_nomag4_mea.Fmag(nomag4_endi),10,'filled','r');

subplot(2,4,5); grid on; hold on;
xlabel('Insertion Depth [mm]'); title('Trial 1');
ylabel('Guided ||Force|| [mN]');
ylim([0 maxY]);
plot(data_mag1.depth_insertion, data_mag1.Fmag, 'Color', [colorsMat(1,:), 0.3*alpha], 'LineWidth', 1);
plot(data_mag1.depth_insertion, data_mag1.Fmag_smooth, 'Color',colorsMat(1,:), 'LineWidth', 2);
scatter(data_mag1.depth_insertion(mag1_endi),data_mag1.Fmag(mag1_endi),10,'filled','r');

subplot(2,4,6); grid on; hold on;
xlabel('Insertion Depth [mm]'); title('Trial 2');
ylim([0 maxY]);
plot(data_mag2.depth_insertion, data_mag2.Fmag, 'Color', [colorsMat(1,:), 0.3*alpha], 'LineWidth', 1);
plot(data_mag2.depth_insertion, data_mag2.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth', 2);
scatter(data_mag2.depth_insertion(mag2_endi),data_mag2.Fmag(mag2_endi),10,'filled','r');

subplot(2,4,7); grid on; hold on;
xlabel('Insertion Depth [mm]'); title('Trial 3');
ylim([0 maxY]);
plot(data_mag3.depth_insertion, data_mag3.Fmag, 'Color', [colorsMat(1,:), 0.3*alpha], 'LineWidth', 1);
plot(data_mag3.depth_insertion, data_mag3.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth', 2);
scatter(data_mag3.depth_insertion(mag3_endi),data_mag3.Fmag(mag3_endi),10,'filled','r');

subplot(2,4,8); grid on; hold on;
xlabel('Insertion Depth [mm]'); title('Trial 3');
ylim([0 maxY]);
plot(data_mag4.depth_insertion, data_mag4.Fmag, 'Color', [colorsMat(1,:), 0.3*alpha], 'LineWidth', 1);
plot(data_mag4.depth_insertion, data_mag4.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth', 2);
scatter(data_mag4.depth_insertion(mag4_endi),data_mag4.Fmag(mag4_endi),10,'filled','r');

%% Plot 2: Averages
maxX = min([data_nomag1_mea.depth_insertion(nomag1_endi);...
            data_nomag2_mea.depth_insertion(nomag2_endi);...
            data_nomag3_mea.depth_insertion(nomag3_endi);...
            data_nomag4_mea.depth_insertion(nomag4_endi);...
            data_mag1.depth_insertion(mag1_endi);...
            data_mag2.depth_insertion(mag2_endi);...
            data_mag3.depth_insertion(mag3_endi);...
            data_mag4.depth_insertion(mag4_endi)]);

xvec = linspace(0.02,maxX,1000);

Fnomag_mea = [interp1(data_nomag1_mea.depth_insertion(1:nomag1_endi),...
                      data_nomag1_mea.Fmag_smooth(1:nomag1_endi),xvec);...
              interp1(data_nomag2_mea.depth_insertion(1:nomag2_endi),...
                      data_nomag2_mea.Fmag_smooth(1:nomag2_endi),xvec);...
              interp1(data_nomag3_mea.depth_insertion(1:nomag3_endi),...
                      data_nomag3_mea.Fmag_smooth(1:nomag3_endi),xvec);...
              interp1(data_nomag4_mea.depth_insertion(1:nomag4_endi),...
                      data_nomag4_mea.Fmag_smooth(1:nomag4_endi),xvec)];
          
Fmag = [interp1(data_mag1.depth_insertion(1:mag1_endi),...
                data_mag1.Fmag_smooth(1:mag1_endi),xvec);...
        interp1(data_mag2.depth_insertion(1:mag2_endi),...
                data_mag2.Fmag_smooth(1:mag2_endi),xvec);...
        interp1(data_mag3.depth_insertion(1:mag3_endi),...
                data_mag3.Fmag_smooth(1:mag3_endi),xvec);...
        interp1(data_mag4.depth_insertion(1:mag4_endi),...
                data_mag4.Fmag_smooth(1:mag4_endi),xvec)];
    
% Fnomag_ea = [interp1(data_nomag1_ea.depth_insertion(2:end-1),...
%                      data_nomag1_ea.Fmag_smooth(2:end-1),xvec);...
%              interp1(data_nomag2_ea.depth_insertion(2:end-1),...
%                      data_nomag2_ea.Fmag_smooth(2:end-1),xvec);...
%              interp1(data_nomag3_ea.depth_insertion(2:end-1),...
%                      data_nomag3_ea.Fmag_smooth(2:end-1),xvec);...
%              interp1(data_nomag4_ea.depth_insertion(2:end-1),...
%                      data_nomag4_ea.Fmag_smooth(2:end-1),xvec)];

% Compute Averages and Standard Deviations
Favg_nomag_mea = mean(Fnomag_mea,1);
std_nomag_mea = std(Fnomag_mea);

Favg_mag = mean(Fmag,1);
std_mag = std(Fmag);

% Favg_nomag_ea = nanmean(Fnomag_ea,1);
% std_nomag_ea = std(Fnomag_ea);

% Plot the averages with no shifting
figure(3); clf(3); grid on; hold on;
xlim([0,xvec(end)]); 
ylim([0,120]);
xlabel('Insertion Depth [mm]'); ylabel('Average ||Force|| [mN]');
h1 = plot(xvec, Favg_nomag_mea, 'Color', 'r',...
    'LineWidth',1,'LineStyle','--');
fill([xvec fliplr(xvec)],[Favg_nomag_mea + std_nomag_mea,...
    fliplr(Favg_nomag_mea-std_nomag_mea)],'r','FaceAlpha',0.2,...
    'LineStyle','none');
h2 = plot(xvec, Favg_mag, 'Color', 'b',...
    'LineWidth',2,'LineStyle','-');
fill([xvec fliplr(xvec)],[Favg_mag + std_mag,...
    fliplr(Favg_mag-std_mag)],'b','FaceAlpha',0.2,'LineStyle','none');

legend([h1,h2],'No Magnetic Guidance N=4',...
    'Magnetic Guidance N=4');

%% Save Raw Force Magnitudes
p_ug_Fmag = [data_nomag1_mea.Fmag(1:nomag1_endi);...
             data_nomag2_mea.Fmag(1:nomag2_endi);...
             data_nomag3_mea.Fmag(1:nomag3_endi);
             data_nomag4_mea.Fmag(1:nomag4_endi)];
p_g_Fmag = [data_mag1.Fmag(1:mag1_endi);...
            data_mag2.Fmag(1:mag2_endi);...
            data_mag3.Fmag(1:mag3_endi);...
            data_mag4.Fmag(1:mag4_endi)];

save('data\phantom\p_ug_Fmag.mat','p_ug_Fmag');
save('data\phantom\p_g_Fmag.mat','p_g_Fmag');

%% Plot XYZ data
if plotXYZData
    figure(4); clf(4);
    sgtitle('XYZ Data from Phantom Insertions');
    %--------
    % Unguided EA
    %---------
    subplot(3,4,1); grid on; hold on;
    ylabel('Unguided XYZ Forces EA [mN]');
    xlim([0 28]); %ylim([-120,40]);
    title('Trial 1');
    plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
    h1 = plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
    plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
    h2 = plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
    plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
    h3 = plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
    legend([h1,h2,h3],'Fx','Fy','Fz');

    subplot(3,4,2); grid on; hold on;
    xlim([0 28]); %ylim([-120,40]);
    title('Trial 2');
    plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
    plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
    plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
    plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
    plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
    plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
    legend([h1,h2,h3],'Fx','Fy','Fz');

    subplot(3,4,3); grid on; hold on;
    xlim([0 28]); %ylim([-120,40]);
    title('Trial 3');
    plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
    plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
    plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
    plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
    plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
    plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
    legend([h1,h2,h3],'Fx','Fy','Fz');

    subplot(3,4,4); grid on; hold on;
    xlim([0 28]); %ylim([-120,40]);
    title('Trial 4');
    plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
    plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
    plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
    plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
    plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
    plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
    legend([h1,h2,h3],'Fx','Fy','Fz');

    %--------
    % Unguided MEA
    %---------

    subplot(3,4,5); grid on; hold on;
    ylabel('Unguided XYZ Forces MEA [mN]');
    xlim([0 28]); %ylim([-120,40]);
    plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
    h1 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
    plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
    h2 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
    plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
    h3 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
    legend([h1,h2,h3],'Fx','Fy','Fz');

    subplot(3,4,6); grid on; hold on;
    xlim([0 28]); %ylim([-120,40]);
    plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
    plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
    plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
    plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
    plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
    plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
    legend([h1,h2,h3],'Fx','Fy','Fz');

    subplot(3,4,7); grid on; hold on;
    xlim([0 28]); %ylim([-120,40]);
    plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
    plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
    plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
    plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
    plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
    plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
    legend([h1,h2,h3],'Fx','Fy','Fz');

    subplot(3,4,8); grid on; hold on;
    xlim([0 28]); %ylim([-120,40]);
    plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
    plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
    plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
    plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
    plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
    plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
    legend([h1,h2,h3],'Fx','Fy','Fz');

    %--------
    % Guided
    %---------

    subplot(3,4,9); grid on; hold on;
    xlim([0 28]); %ylim([-120,40]);
    xlabel('Insertion Depth (mm)'); ylabel('Guided XYZ Forces [mN]');
    plot(data_mag1.depth_insertion, data_mag1.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
    plot(data_mag1.depth_insertion, data_mag1.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
    plot(data_mag1.depth_insertion, data_mag1.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
    plot(data_mag1.depth_insertion, data_mag1.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
    plot(data_mag1.depth_insertion, data_mag1.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
    plot(data_mag1.depth_insertion, data_mag1.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
    legend([h1,h2,h3],'Fx','Fy','Fz');

    subplot(3,4,10); grid on; hold on;
    xlim([0 28]); %ylim([-120,40]);
    xlabel('Insertion Depth (mm)');
    plot(data_mag2.depth_insertion, data_mag2.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
    plot(data_mag2.depth_insertion, data_mag2.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
    plot(data_mag2.depth_insertion, data_mag2.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
    plot(data_mag2.depth_insertion, data_mag2.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
    plot(data_mag2.depth_insertion, data_mag2.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
    plot(data_mag2.depth_insertion, data_mag2.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
    legend([h1,h2,h3],'Fx','Fy','Fz');

    subplot(3,4,11); grid on; hold on;
    xlim([0 28]); %ylim([-120,40]);
    xlabel('Insertion Depth (mm)');
    plot(data_mag3.depth_insertion, data_mag3.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
    plot(data_mag3.depth_insertion, data_mag3.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
    plot(data_mag3.depth_insertion, data_mag3.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
    plot(data_mag3.depth_insertion, data_mag3.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
    plot(data_mag3.depth_insertion, data_mag3.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
    plot(data_mag3.depth_insertion, data_mag3.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
    legend([h1,h2,h3],'Fx','Fy','Fz');

    subplot(3,4,12); grid on; hold on;
    xlim([0 28]); %ylim([-120,40]);
    xlabel('Insertion Depth (mm)');
    plot(data_mag4.depth_insertion, data_mag4.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
    plot(data_mag4.depth_insertion, data_mag4.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
    plot(data_mag4.depth_insertion, data_mag4.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
    plot(data_mag4.depth_insertion, data_mag4.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
    plot(data_mag4.depth_insertion, data_mag4.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
    plot(data_mag4.depth_insertion, data_mag4.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
    legend([h1,h2,h3],'Fx','Fy','Fz');
end
