%% import force and smaract data
clear all; close all; clc;

addpath('classDef','functions','data');

% filepaths to CSVs exported from ROS bags
base_path = 'data\phantom\';

filepaths_nomag1.smaract = fullfile(base_path, 'UG-MEA\phantom_ug_mea1_trial1_smaract.csv');
filepaths_nomag1.force   = fullfile(base_path, 'UG-MEA\phantom_ug_mea1_trial1_force.csv');

filepaths_nomag2.smaract = fullfile(base_path, 'UG-MEA\phantom_ug_mea1_trial2_smaract.csv');
filepaths_nomag2.force   = fullfile(base_path, 'UG-MEA\phantom_ug_mea1_trial2_force.csv');

filepaths_nomag3.smaract = fullfile(base_path, 'UG-MEA\phantom_ug_mea1_trial3_smaract.csv');
filepaths_nomag3.force   = fullfile(base_path, 'UG-MEA\phantom_ug_mea1_trial3_force.csv');

filepaths_mag1.smaract   = fullfile(base_path, 'guided\phantom_g_mea1_trial1_smaract.csv');
filepaths_mag1.force     = fullfile(base_path, 'guided\phantom_g_mea1_trial1_force.csv');

filepaths_mag2.smaract = fullfile(base_path, 'guided\phantom_g_mea1_trial2_smaract.csv');
filepaths_mag2.force   = fullfile(base_path, 'guided\phantom_g_mea1_trial2_force.csv');

filepaths_mag3.smaract = fullfile(base_path, 'guided\phantom_g_flex31_trimmed_smaract.csv');
filepaths_mag3.force   = fullfile(base_path, 'guided\phantom_g_flex31_trimmed_force.csv');

filepaths_mag4.smaract = fullfile(base_path, 'guided\phantom_g_mea1_trial3_trimmed_smaract.csv');
filepaths_mag4.force   = fullfile(base_path, 'guided\phantom_g_mea1_trial3_trimmed_force.csv');

filepaths_mag5.smaract = fullfile(base_path, 'guided\phantom_g_mea1_trial4_smaract.csv');
filepaths_mag5.force   = fullfile(base_path, 'guided\phantom_g_mea1_trial4_force.csv');

filepaths_mag6.smaract = fullfile(base_path, 'guided\phantom_g_mea1_trial5_smaract.csv');
filepaths_mag6.force   = fullfile(base_path, 'guided\phantom_g_mea1_trial5_force.csv');

filepaths_mag7.smaract = fullfile(base_path, 'guided\phantom_g_mea1_trial6_smaract.csv');
filepaths_mag7.force   = fullfile(base_path, 'guided\phantom_g_mea1_trial6_force.csv');

filepaths_mag8.smaract = fullfile(base_path, 'guided\phantom_g_mea1_trial7_smaract.csv');
filepaths_mag8.force   = fullfile(base_path, 'guided\phantom_g_mea1_trial7_force.csv');

% create data objects
data_nomag1 = MagneticGuidanceData(filepaths_nomag1);
data_nomag2 = MagneticGuidanceData(filepaths_nomag2);
data_nomag3 = MagneticGuidanceData(filepaths_nomag3);

data_mag1 = MagneticGuidanceData(filepaths_mag1);
data_mag2 = MagneticGuidanceData(filepaths_mag2);
data_mag3 = MagneticGuidanceData(filepaths_mag3);
data_mag4 = MagneticGuidanceData(filepaths_mag4);
data_mag5 = MagneticGuidanceData(filepaths_mag5);
data_mag6 = MagneticGuidanceData(filepaths_mag6);
data_mag7 = MagneticGuidanceData(filepaths_mag7);
data_mag8 = MagneticGuidanceData(filepaths_mag8);

% set smoothing span (proportion of data points)
beta = 0.01;

data_nomag1.setSmoothSpan(beta);
data_nomag2.setSmoothSpan(beta);
data_nomag3.setSmoothSpan(beta);

data_mag1.setSmoothSpan(beta);
data_mag2.setSmoothSpan(beta);
data_mag3.setSmoothSpan(beta);
data_mag4.setSmoothSpan(beta);
data_mag5.setSmoothSpan(beta);
data_mag6.setSmoothSpan(beta);
data_mag7.setSmoothSpan(beta);
data_mag8.setSmoothSpan(beta);

%% plot
PlotColors; % load custom colors
alpha = 1; % reduce transparency of unguided plot lines

% Figure: force components
figure(1); clf(1); hold on; grid on;
% set(gca,'FontSize',15)
% xlim([data_nomag.depth_insertion(1)-0.5, data_nomag.depth_insertion(end)+0.5])
title('Guided vs Unguided Insertion')
xlabel('Insertion Depth [mm]')
ylabel('Force [mN]')

colors1 = [colors.green;colors.red;colors.purple];
colors2 = ['k';'b';'m'];
colors3 = [colors.brown;colors.light_orange;colors.light_green];
colors4 = [colors.gray;colors.orange;colors.light_blue];
colors5 = [colors.light_red;colors.orange;colors.brown];
    
h1 = addData2Plot(data_nomag1,alpha,colors1);
% h2 = addData2Plot(data_mag1,alpha,colors2);
% h3 = addData2Plot(data_nomag2,alpha,colors3);
% h4 = addData2Plot(data_mag2,alpha,colors4);
h5 = addData2Plot(data_mag3,alpha,colors2);
% h6 = addData2Plot(data_mag4,alpha,colors3);

[~,icons,~,~] = legend([h1(1), h5(1), h1(2), h5(2), h1(3), h5(3)], ...
        'Fx_{G}Flex28','Fx_{G}Flex31','Fy_{G}Flex28','Fy_{G}Flex31',...
        'Fz_{G}Flex28','Fz_{G}Flex31', 'Location','southwest');

% Figure: force magnitude
figure(2); clf(2); hold on; grid on;

% plot(data_nomag1.depth_insertion, data_nomag1.Fmag, 'Color', [colors1(1,:),  0.3*alpha], 'LineWidth',1);
h1 = plot(data_nomag1.depth_insertion, data_nomag1.Fmag_smooth, 'Color', [colors1(1,:),  alpha], 'LineWidth',2);

% plot(data_nomag2.depth_insertion, data_nomag2.Fmag, 'Color', [colors3(1,:),  0.3*alpha], 'LineWidth',1);
h2 = plot(data_nomag2.depth_insertion, data_nomag2.Fmag_smooth, 'Color', [colors3(1,:),  alpha], 'LineWidth',2);

% plot(data_nomag2.depth_insertion, data_nomag2.Fmag, 'Color', [colors3(1,:),  0.3*alpha], 'LineWidth',1);
h3 = plot(data_nomag3.depth_insertion, data_nomag3.Fmag_smooth, 'Color', 'b', 'LineWidth',2);

% plot(data_mag1.depth_insertion, data_mag1.Fmag, 'Color', [colors2(1,:), 0.3*alpha], 'LineWidth', 1);
h4 = plot(data_mag1.depth_insertion, data_mag1.Fmag_smooth, 'Color', colors2(1,:), 'LineWidth', 2);

% plot(data_mag2.depth_insertion, data_mag2.Fmag, 'Color', [colors4(1,:), 0.3*alpha], 'LineWidth', 1);
h5 = plot(data_mag2.depth_insertion, data_mag2.Fmag_smooth, 'Color', colors4(1,:), 'LineWidth', 2);

% plot(data_mag3.depth_insertion, data_mag3.Fmag, 'Color', [colors5(1,:), 0.3*alpha], 'LineWidth', 1);
h6 = plot(data_mag3.depth_insertion, data_mag3.Fmag_smooth, 'Color', colors5(1,:), 'LineWidth', 2);

% plot(data_mag4.depth_insertion, data_mag4.Fmag, 'Color', ['g', 0.3*alpha], 'LineWidth', 1);
h7 = plot(data_mag4.depth_insertion, data_mag4.Fmag_smooth, 'Color', 'g', 'LineWidth', 2);

% plot(data_mag5.depth_insertion, data_mag5.Fmag, 'Color', ['r', 0.3*alpha], 'LineWidth', 1);
h8 = plot(data_mag5.depth_insertion, data_mag5.Fmag_smooth, 'Color', 'r', 'LineWidth', 2);

% plot(data_mag6.depth_insertion, data_mag6.Fmag, 'Color', ['y', 0.3*alpha], 'LineWidth', 1);
h9 = plot(data_mag6.depth_insertion, data_mag6.Fmag_smooth, 'Color', 'y', 'LineWidth', 2);

% plot(data_mag7.depth_insertion, data_mag7.Fmag, 'Color', ['m', 0.3*alpha], 'LineWidth', 1);
h10 = plot(data_mag7.depth_insertion, data_mag7.Fmag_smooth, 'Color', 'm', 'LineWidth', 2);

% plot(data_mag8.depth_insertion, data_mag8.Fmag, 'Color', ['m', 0.3*alpha], 'LineWidth', 1);
h11 = plot(data_mag8.depth_insertion+1, data_mag8.Fmag_smooth, 'Color', rand(1,3), 'LineWidth', 2);

% nomag_Fmag_interp = interp1(data_nomag1.depth_insertion, data_nomag1.Fmag_smooth, data_mag.depth_insertion);
% plot(data_mag.depth_insertion, (data_mag.Fmag_smooth - nomag_Fmag_interp), 'Color', colors.blue, 'LineWidth', 2);

% xlim([data_nomag.depth_insertion(1)-0.5, data_nomag.depth_insertion(end)+0.5])
title('Guided vs Unguided Insertion')
xlabel('Insertion Depth [mm]')
ylabel('||Force|| [mN]')
legend([h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11],{'nomag1','nomag2','nomag3','mag1','mag2','mag3(flex31)','mag4','mag5','mag6','mag7','mag8'});
xlim([0,28]);
ylim([0,40]);