%% Data Analysis for Phantom Trials
% Trevor Bruns and Katy Riojas
% Last Updated: 9/5/19

clear all; clc;

%% INPUTS
addpath('classDef','functions','data');
load('avg_cal_slope.mat');
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

% Create data objects
data_nomag1_mea = MagneticGuidanceData(filepaths_nomag1_mea,avg_cal_slope);
data_nomag2_mea = MagneticGuidanceData(filepaths_nomag2_mea,avg_cal_slope);
data_nomag3_mea = MagneticGuidanceData(filepaths_nomag3_mea,avg_cal_slope);
data_nomag4_mea = MagneticGuidanceData(filepaths_nomag4_mea,avg_cal_slope);

data_mag1 = MagneticGuidanceData(filepaths_mag1,avg_cal_slope);
data_mag2 = MagneticGuidanceData(filepaths_mag2,avg_cal_slope);
data_mag3 = MagneticGuidanceData(filepaths_mag3,avg_cal_slope);
data_mag4 = MagneticGuidanceData(filepaths_mag4,avg_cal_slope);

data_nomag1_ea = MagneticGuidanceData(filepaths_nomag1_ea,avg_cal_slope);
data_nomag2_ea = MagneticGuidanceData(filepaths_nomag2_ea,avg_cal_slope);
data_nomag3_ea = MagneticGuidanceData(filepaths_nomag3_ea,avg_cal_slope);
data_nomag4_ea = MagneticGuidanceData(filepaths_nomag4_ea,avg_cal_slope);

% Smooth Data
data_nomag1_ea = data_nomag1_ea.setSmoothSpan(smooth_span);
data_nomag2_ea = data_nomag2_ea.setSmoothSpan(smooth_span);
data_nomag3_ea = data_nomag3_ea.setSmoothSpan(smooth_span);
data_nomag4_ea = data_nomag4_ea.setSmoothSpan(smooth_span);

data_nomag1_mea = data_nomag1_mea.setSmoothSpan(smooth_span);
data_nomag2_mea = data_nomag2_mea.setSmoothSpan(smooth_span);
data_nomag3_mea = data_nomag3_mea.setSmoothSpan(smooth_span);
data_nomag4_mea = data_nomag4_mea.setSmoothSpan(smooth_span);

data_mag1 = data_mag1.setSmoothSpan(smooth_span);
data_mag2 = data_mag2.setSmoothSpan(smooth_span);
data_mag3 = data_mag3.setSmoothSpan(smooth_span);
data_mag4 = data_mag4.setSmoothSpan(smooth_span);

% Calculate Rotated Forces
data_nomag1_ea = data_nomag1_ea.setT(T_fixture_st);
data_nomag2_ea = data_nomag2_ea.setT(T_fixture_st);
data_nomag3_ea = data_nomag3_ea.setT(T_fixture_st);
data_nomag4_ea = data_nomag4_ea.setT(T_fixture_st);

data_nomag1_mea = data_nomag1_mea.setT(T_fixture_st);
data_nomag2_mea = data_nomag2_mea.setT(T_fixture_st);
data_nomag3_mea = data_nomag3_mea.setT(T_fixture_st);
data_nomag4_mea = data_nomag4_mea.setT(T_fixture_st);

data_mag1 = data_mag1.setT(T_fixture_st);
data_mag2 = data_mag2.setT(T_fixture_st);
data_mag3 = data_mag3.setT(T_fixture_st);
data_mag4 = data_mag4.setT(T_fixture_st);

%% Plot 1: XYZ forces (smoothed and raw) vs. insertion depth
alpha = 1; % line opacity
colorsMat = distinguishable_colors(12);
xyzColor = distinguishable_colors(3); % xyz colors

figure(1); clf(1);
sgtitle('Phantom Insertions');
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

%% Figure 2: Force Magnitude Trial Comparison
figure(2); clf(2); hold on; grid on;
title('Guided vs Unguided Insertion');
xlabel('Insertion Depth [mm]');
ylabel('||Force|| [mN]');
% --------------
% Unguided EA
% --------------
plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fmag, 'Color', [colorsMat(1,:),  0.3*alpha], 'LineWidth',1);
h1 = plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth',3,'LineStyle',':');
 
plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fmag, 'Color', colorsMat(2,:), 'LineWidth',1);
h2 = plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fmag_smooth, 'Color', colorsMat(2,:), 'LineWidth',3,'LineStyle',':');

plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fmag, 'Color', [colorsMat(3,:),  0.3*alpha], 'LineWidth',1);
h3 = plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fmag_smooth, 'Color', colorsMat(3,:), 'LineWidth',3,'LineStyle',':');

plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fmag, 'Color', [colorsMat(4,:),  0.3*alpha], 'LineWidth',1);
h4 = plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fmag_smooth, 'Color', colorsMat(4,:), 'LineWidth',3,'LineStyle',':');

% --------------
% Unguided MEA
% --------------

plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fmag, 'Color', [colorsMat(5,:),  0.3*alpha], 'LineWidth',1);
h5 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fmag_smooth, 'Color', colorsMat(5,:), 'LineWidth',1,'LineStyle','--');

plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fmag, 'Color', colorsMat(6,:), 'LineWidth',1);
h6 = plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fmag_smooth, 'Color', colorsMat(6,:), 'LineWidth',1,'LineStyle','--');

plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fmag, 'Color', [colorsMat(7,:),  0.3*alpha], 'LineWidth',1);
h7 = plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fmag_smooth, 'Color', colorsMat(7,:), 'LineWidth',1,'LineStyle','--');

plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fmag, 'Color', [colorsMat(8,:),  0.3*alpha], 'LineWidth',1);
h8 = plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fmag_smooth, 'Color', colorsMat(8,:), 'LineWidth',1,'LineStyle','--');

% --------------
% Guided
% --------------
plot(data_mag1.depth_insertion, data_mag1.Fmag, 'Color', [colorsMat(9,:), 0.3*alpha], 'LineWidth', 1);
h9 = plot(data_mag1.depth_insertion, data_mag1.Fmag_smooth, 'Color',colorsMat(9,:), 'LineWidth', 2);

plot(data_mag2.depth_insertion, data_mag2.Fmag, 'Color', [colorsMat(10,:), 0.3*alpha], 'LineWidth', 1);
h10 = plot(data_mag2.depth_insertion, data_mag2.Fmag_smooth, 'Color', colorsMat(10,:), 'LineWidth', 2);

plot(data_mag3.depth_insertion, data_mag3.Fmag, 'Color', [colorsMat(11,:), 0.3*alpha], 'LineWidth', 1);
h11 = plot(data_mag3.depth_insertion, data_mag3.Fmag_smooth, 'Color', colorsMat(11,:), 'LineWidth', 2);

plot(data_mag4.depth_insertion, data_mag4.Fmag, 'Color', [colorsMat(12,:), 0.3*alpha], 'LineWidth', 1);
h12 = plot(data_mag4.depth_insertion, data_mag4.Fmag_smooth, 'Color', colorsMat(12,:), 'LineWidth', 2);

legend([h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12],{'nomag1-ea',...
    'nomag2-ea','nomag3-ea','nomag4-ea','nomag1-mea',...
    'nomag2-mea','nomag3-mea','nomag4-mea','mag1','mag2','mag3','mag4'});

%% Plot 3: Averages
xvec = linspace(0.25,23.7,1000);

Fnomag_mea = [interp1(data_nomag1_mea.depth_insertion(2:end-1),...
                      data_nomag1_mea.Fmag_smooth(2:end-1),xvec);...
              interp1(data_nomag2_mea.depth_insertion(2:end-1),...
                      data_nomag2_mea.Fmag_smooth(2:end-1),xvec);...
              interp1(data_nomag3_mea.depth_insertion(2:end-1),...
                      data_nomag3_mea.Fmag_smooth(2:end-1),xvec);...
              interp1(data_nomag4_mea.depth_insertion(2:end-1),...
                      data_nomag4_mea.Fmag_smooth(2:end-1),xvec)];
          
Fmag = [interp1(data_mag1.depth_insertion(2:end-1),...
                data_mag1.Fmag_smooth(2:end-1),xvec);...
        interp1(data_mag2.depth_insertion(2:end-1),...
                data_mag2.Fmag_smooth(2:end-1),xvec);...
        interp1(data_mag3.depth_insertion(2:end-1),...
                data_mag3.Fmag_smooth(2:end-1),xvec);...
        interp1(data_mag4.depth_insertion(2:end-1),...
                data_mag4.Fmag_smooth(2:end-1),xvec)];
    
Fnomag_ea = [interp1(data_nomag1_ea.depth_insertion(2:end-1),...
                     data_nomag1_ea.Fmag_smooth(2:end-1),xvec);...
             interp1(data_nomag2_ea.depth_insertion(2:end-1),...
                     data_nomag2_ea.Fmag_smooth(2:end-1),xvec);...
             interp1(data_nomag3_ea.depth_insertion(2:end-1),...
                     data_nomag3_ea.Fmag_smooth(2:end-1),xvec);...
             interp1(data_nomag4_ea.depth_insertion(2:end-1),...
                     data_nomag4_ea.Fmag_smooth(2:end-1),xvec)];

% Compute Averages and Standard Deviations
Favg_nomag_mea = mean(Fnomag_mea,1);
std_nomag_mea = std(Fnomag_mea);

Favg_mag = mean(Fmag,1);
std_mag = std(Fmag);

Favg_nomag_ea = nanmean(Fnomag_ea,1);
std_nomag_ea = std(Fnomag_ea);

% Plot the averages with no shifting
figure(3); clf(3); grid on; hold on;
xlim([0,xvec(end)]); ylim([0,120]);
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
% h3 = plot(xvec,-(Favg_nomag_mea-Favg_mag),'Color','k',...
%     'LineWidth',4,'LineStyle',':');

legend([h1,h2],'No Magnetic Guidance N=4',...
    'Magnetic Guidance N=4');%,'Average Difference N=4');


%% Plot 4: XYZ Forces in ST frame
figure(4); clf(4);
sgtitle('Phantom Insertions XYZ Forces');

%--------
subplot(3,4,1); grid on; hold on;
ylabel('Unguided XYZ Forces EA [mN]');
xlim([0 28]); ylim([-75,250]);
title('Trial 1');
% plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fx_smooth_st, 'Color', [xyzColor(1,:), 0.3*alpha]);
h1 = plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
h2 = plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
h3 = plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);
% legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(3,4,2); grid on; hold on;
xlim([0 28]); ylim([-75,250]);
title('Trial 2');
% plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);
% legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(3,4,3); grid on; hold on;
xlim([0 28]); ylim([-75,250]);
title('Trial 3');
% plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);
% legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(3,4,4); grid on; hold on;
xlim([0 28]); ylim([-75,250]);
title('Trial 4');
% plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);
% legend([h1,h2,h3],'Fx','Fy','Fz');

%--------

subplot(3,4,5); grid on; hold on;
ylabel('Unguided XYZ Forces MEA [mN]');
xlim([0 28]); ylim([-75,250]);
% plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
h1 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
h2 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
h3 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);
% legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(3,4,6); grid on; hold on;
xlim([0 28]); ylim([-75,250]);
% plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);
% legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(3,4,7); grid on; hold on;
xlim([0 28]); ylim([-75,250]);
% plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);
% legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(3,4,8); grid on; hold on;
xlim([0 28]); ylim([-75,250]);
% plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);
% legend([h1,h2,h3],'Fx','Fy','Fz');

%--------

subplot(3,4,9); grid on; hold on;
xlim([0 28]); ylim([-75,250]);
xlabel('Insertion Depth (mm)'); ylabel('Guided XYZ Forces [mN]');
% plot(data_mag1.depth_insertion, data_mag1.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_mag1.depth_insertion, data_mag1.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_mag1.depth_insertion, data_mag1.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_mag1.depth_insertion, data_mag1.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_mag1.depth_insertion, data_mag1.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_mag1.depth_insertion, data_mag1.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);
% legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(3,4,10); grid on; hold on;
xlim([0 28]); ylim([-75,250]);
xlabel('Insertion Depth (mm)');
% plot(data_mag2.depth_insertion, data_mag2.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_mag2.depth_insertion, data_mag2.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_mag2.depth_insertion, data_mag2.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_mag2.depth_insertion, data_mag2.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_mag2.depth_insertion, data_mag2.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_mag2.depth_insertion, data_mag2.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);
% legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(3,4,11); grid on; hold on;
xlim([0 28]); ylim([-75,250]);
xlabel('Insertion Depth (mm)');
% plot(data_mag3.depth_insertion, data_mag3.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_mag3.depth_insertion, data_mag3.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_mag3.depth_insertion, data_mag3.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_mag3.depth_insertion, data_mag3.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_mag3.depth_insertion, data_mag3.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_mag3.depth_insertion, data_mag3.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);
% legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(3,4,12); grid on; hold on;
xlim([0 28]); ylim([-75,250]);
xlabel('Insertion Depth (mm)');
% plot(data_mag4.depth_insertion, data_mag4.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_mag4.depth_insertion, data_mag4.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_mag4.depth_insertion, data_mag4.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_mag4.depth_insertion, data_mag4.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_mag4.depth_insertion, data_mag4.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_mag4.depth_insertion, data_mag4.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);
legend([h1,h2,h3],'-Fz','Fx','-Fy');

%% Save Raw Force Magnitudes
p_ug_Fmag = [data_nomag1_mea.Fmag;data_nomag2_mea.Fmag;...
             data_nomag3_mea.Fmag;data_nomag4_mea.Fmag];
p_g_Fmag = [data_mag1.Fmag;data_mag2.Fmag;data_mag3.Fmag;data_mag4.Fmag];

save('data\phantom\p_ug_Fmag.mat','p_ug_Fmag');
save('data\phantom\p_g_Fmag.mat','p_g_Fmag');
