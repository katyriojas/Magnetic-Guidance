%% import force and smaract data
clear all; close all; clc;

addpath('classDef','functions','data');

% pull in the 1.25mm/s guided mea2 data

% filepaths to CSVs exported from ROS bags
base_path = 'C:\Users\riojaske\Documents\magsteer\Magnetic-Guidance\data\phantom\';

filepaths_nomag1.smaract = fullfile(base_path, 'unguided_saline_1.25\phantom_ug_mea1_trial1_1.25_smaract.csv');
filepaths_nomag1.force   = fullfile(base_path, 'unguided_saline_1.25\phantom_ug_mea1_trial1_1.25_force.csv');

filepaths_nomag2.smaract = fullfile(base_path, 'unguided_saline_1.25\phantom_ug_mea1_trial2_1.25_smaract.csv');
filepaths_nomag2.force   = fullfile(base_path, 'unguided_saline_1.25\phantom_ug_mea1_trial2_1.25_force.csv');

filepaths_nomag3.smaract = fullfile(base_path, 'unguided_saline_1.25\phantom_ug_mea1_trial3_1.25_smaract.csv');
filepaths_nomag3.force   = fullfile(base_path, 'unguided_saline_1.25\phantom_ug_mea1_trial3_1.25_force.csv');

filepaths_nomag4.smaract = fullfile(base_path, 'unguided_saline_1.25\phantom_ug_mea1_trial4_1.25_smaract.csv');
filepaths_nomag4.force   = fullfile(base_path, 'unguided_saline_1.25\phantom_ug_mea1_trial4_1.25_force.csv');

filepaths_mag1.smaract   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial1_1.25_smaract.csv');
filepaths_mag1.force     = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial1_1.25_force.csv');

filepaths_mag2.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial2_1.25_smaract.csv');
filepaths_mag2.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial2_1.25_force.csv');

filepaths_mag3.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial3_1.25_smaract.csv');
filepaths_mag3.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial3_1.25_force.csv');

filepaths_mag4.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial4_1.25_smaract.csv');
filepaths_mag4.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial4_1.25_force.csv');

% create data objects
data_nomag1 = MagneticGuidanceData(filepaths_nomag1);
data_nomag2 = MagneticGuidanceData(filepaths_nomag2);
data_nomag3 = MagneticGuidanceData(filepaths_nomag3);
data_nomag4 = MagneticGuidanceData(filepaths_nomag4);

data_mag1 = MagneticGuidanceData(filepaths_mag1);
data_mag2 = MagneticGuidanceData(filepaths_mag2);
data_mag3 = MagneticGuidanceData(filepaths_mag3);
data_mag4 = MagneticGuidanceData(filepaths_mag4);

% set smoothing span (proportion of data points)
beta = 0.01;

data_nomag1.setSmoothSpan(beta);
data_nomag2.setSmoothSpan(beta);
data_nomag3.setSmoothSpan(beta);
data_nomag4.setSmoothSpan(beta);

data_mag1.setSmoothSpan(beta);
data_mag2.setSmoothSpan(beta);
data_mag3.setSmoothSpan(beta);
data_mag4.setSmoothSpan(beta);

%% plot
PlotColors; % load custom colors
alpha = 1; % reduce transparency of unguided plot lines
% colors1 = [colors.red;colors.red;colors.purple];
% colors2 = ['k';'b';'m'];
% colors3 = [colors.green;colors.light_orange;colors.light_green];
% colors4 = [colors.gray;colors.orange;colors.light_blue];
% colors5 = [colors.red;colors.orange;colors.brown];
% colors6 = [colors.orange;colors.red;colors.light_blue];
% colors7 = 'b';
% colors8 = 'm';

colorsMat = distinguishable_colors(11);

% Figure: force magnitude
figure(1); clf(1); hold on; grid on;

% plot(data_nomag1.depth_insertion, data_nomag1.Fmag, 'Color', [rand(1,3),  0.3*alpha], 'LineWidth',1);
h1 = plot(data_nomag1.depth_insertion, data_nomag1.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth',1,'LineStyle','--');
% 
% plot(data_nomag2.depth_insertion, data_nomag2.Fmag, 'Color', colors2(1,:), 'LineWidth',1);
h2 = plot(data_nomag2.depth_insertion, data_nomag2.Fmag_smooth, 'Color', colorsMat(2,:), 'LineWidth',1,'LineStyle','--');

% plot(data_nomag3.depth_insertion, data_nomag3.Fmag, 'Color', [colors3(1,:),  0.3*alpha], 'LineWidth',1);
h3 = plot(data_nomag3.depth_insertion, data_nomag3.Fmag_smooth, 'Color', colorsMat(3,:), 'LineWidth',1,'LineStyle','--');

% plot(data_nomag4.depth_insertion, data_nomag4.Fmag, 'Color', [colors4(1,:),  0.3*alpha], 'LineWidth',1);
h4 = plot(data_nomag4.depth_insertion, data_nomag4.Fmag_smooth, 'Color', colorsMat(4,:), 'LineWidth',1,'LineStyle','--');

% plot(data_mag1.depth_insertion, data_mag1.Fmag, 'Color', [colors5(1,:), 0.3*alpha], 'LineWidth', 1);
h5 = plot(data_mag1.depth_insertion, data_mag1.Fmag_smooth, 'Color',colorsMat(5,:), 'LineWidth', 2);

% plot(data_mag2.depth_insertion, data_mag2.Fmag, 'Color', [colors6(1,:), 0.3*alpha], 'LineWidth', 1);
h6 = plot(data_mag2.depth_insertion, data_mag2.Fmag_smooth, 'Color', colorsMat(6,:), 'LineWidth', 2);

% plot(data_mag3.depth_insertion, data_mag3.Fmag, 'Color', [colors7(1,:), 0.3*alpha], 'LineWidth', 1);
h7 = plot(data_mag3.depth_insertion, data_mag3.Fmag_smooth, 'Color', colorsMat(7,:), 'LineWidth', 2);

% plot(data_mag4.depth_insertion, data_mag4.Fmag, 'Color', [colors8(1,:), 0.3*alpha], 'LineWidth', 1);
h8 = plot(data_mag4.depth_insertion, data_mag4.Fmag_smooth, 'Color', colorsMat(8,:), 'LineWidth', 2);

% nomag_Fmag_interp = interp1(data_nomag1.depth_insertion, data_nomag1.Fmag_smooth, data_mag1.depth_insertion);
% plot(data_mag1.depth_insertion, (data_mag1.Fmag_smooth - nomag_Fmag_interp), 'Color', colors.blue, 'LineWidth', 2);

title('Guided vs Unguided Insertion')
xlabel('Insertion Depth [mm]')
ylabel('||Force|| [mN]')
% legend([h1,h4],{'nomag1','mag1'});
legend([h1,h2,h3,h4,h5,h6,h7,h8],{'nomag1','nomag2','nomag3','nomag4','mag1','mag2','mag3','mag4'});
xlim([0,27]);
ylim([0,40]);

xvec = linspace(0,23,1000);
Fnomag = [interp1(data_nomag1.depth_insertion(2:end-1),data_nomag1.Fmag_smooth(2:end-1),xvec);...
          interp1(data_nomag2.depth_insertion(2:end-1),data_nomag2.Fmag_smooth(2:end-1),xvec);...
          interp1(data_nomag3.depth_insertion(2:end-1),data_nomag3.Fmag_smooth(2:end-1),xvec);...
          interp1(data_nomag4.depth_insertion(2:end-1),data_nomag4.Fmag_smooth(2:end-1),xvec)];
Fmag = [interp1(data_mag1.depth_insertion(2:end-1),data_mag1.Fmag_smooth(2:end-1),xvec);...
        interp1(data_mag2.depth_insertion(2:end-1),data_mag2.Fmag_smooth(2:end-1),xvec);...
        interp1(data_mag3.depth_insertion(2:end-1),data_mag3.Fmag_smooth(2:end-1),xvec);...
        interp1(data_mag4.depth_insertion(2:end-1),data_mag4.Fmag_smooth(2:end-1),xvec)];

Favg_nomag = mean(Fnomag,1);
Favg_mag = mean(Fmag,1);

% figure(2); grid on; hold on;
plot(xvec, Favg_nomag, 'Color', colorsMat(9,:), 'LineWidth',4,'LineStyle',':');
plot(xvec, Favg_mag, 'Color', colorsMat(10,:), 'LineWidth',4,'LineStyle',':');
plot(xvec,-(Favg_nomag-Favg_mag),'Color',colorsMat(11,:),'LineWidth',4,'LineStyle',':');