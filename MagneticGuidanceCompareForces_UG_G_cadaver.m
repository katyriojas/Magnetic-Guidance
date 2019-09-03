%% import force and smaract data
clear all; close all; clc;

addpath('classDef','functions','data');

% Pull in the 1.25mm/s guided mea2 data

% Filepaths to CSVs exported from ROS bags
base_path = 'C:\Users\riojaske\Documents\magsteer\Magnetic-Guidance\data\cadaver';

filepaths_nomag1_mea.smaract = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial1_smaract.csv');
filepaths_nomag1_mea.force   = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial1_force.csv');

filepaths_nomag2_mea.smaract = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial2_smaract.csv');
filepaths_nomag2_mea.force   = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial2_force.csv');
% 
filepaths_nomag3_mea.smaract = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial3_smaract.csv');
filepaths_nomag3_mea.force   = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial3_force.csv');

filepaths_mag1.smaract   = fullfile(base_path, 'guided\cadaver_g_mea2_trial1_smaract.csv');
filepaths_mag1.force     = fullfile(base_path, 'guided\cadaver_g_mea2_trial1_force.csv');

filepaths_mag2.smaract = fullfile(base_path, 'guided\cadaver_g_mea2_trial2_smaract.csv');
filepaths_mag2.force   = fullfile(base_path, 'guided\cadaver_g_mea2_trial2_force.csv');

filepaths_mag3.smaract = fullfile(base_path, 'guided\cadaver_g_mea2_trial3_smaract.csv');
filepaths_mag3.force   = fullfile(base_path, 'guided\cadaver_g_mea2_trial3_force.csv');

% Create data objects
data_nomag1_mea = MagneticGuidanceData(filepaths_nomag1_mea);
data_nomag2_mea = MagneticGuidanceData(filepaths_nomag2_mea);
data_nomag3_mea = MagneticGuidanceData(filepaths_nomag3_mea);

data_mag1 = MagneticGuidanceData(filepaths_mag1);
data_mag2 = MagneticGuidanceData(filepaths_mag2);
data_mag3 = MagneticGuidanceData(filepaths_mag3);

% set smoothing span (proportion of data points)
beta = 0.01;

data_nomag1_mea.setSmoothSpan(beta);
data_nomag2_mea.setSmoothSpan(beta);
data_nomag3_mea.setSmoothSpan(beta);

data_mag1.setSmoothSpan(beta);
data_mag2.setSmoothSpan(beta);
data_mag3.setSmoothSpan(beta);

%% plot

PlotColors; % load custom colors
alpha = 1; % reduce transparency of unguided plot lines

colorsMat = distinguishable_colors(12);

% Figure: force magnitude
figure(1); clf(1); hold on; grid on;

plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fmag, 'Color', [colorsMat(1,:),  0.3*alpha], 'LineWidth',1);
h1 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth',1,'LineStyle','--');

plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fmag, 'Color', [colorsMat(2,:), 0.3*alpha], 'LineWidth',1);
h2 = plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fmag_smooth, 'Color', colorsMat(2,:), 'LineWidth',1,'LineStyle','--');

plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fmag, 'Color', [colorsMat(3,:),  0.3*alpha], 'LineWidth',1);
h3 = plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fmag_smooth, 'Color', colorsMat(3,:), 'LineWidth',1,'LineStyle','--');

plot(data_mag1.depth_insertion, data_mag1.Fmag, 'Color', [colorsMat(4,:), 0.3*alpha], 'LineWidth', 1);
h4 = plot(data_mag1.depth_insertion, data_mag1.Fmag_smooth, 'Color',colorsMat(4,:), 'LineWidth', 2);

plot(data_mag2.depth_insertion, data_mag2.Fmag, 'Color', [colorsMat(5,:), 0.3*alpha], 'LineWidth', 1);
h5 = plot(data_mag2.depth_insertion, data_mag2.Fmag_smooth, 'Color', colorsMat(5,:), 'LineWidth', 2);

plot(data_mag3.depth_insertion, data_mag3.Fmag, 'Color', [colorsMat(6,:), 0.3*alpha], 'LineWidth', 1);
h6 = plot(data_mag3.depth_insertion, data_mag3.Fmag_smooth, 'Color', colorsMat(6,:), 'LineWidth', 2);

title('Guided vs Unguided Insertion')
xlabel('Insertion Depth [mm]')
ylabel('||Force|| [mN]')
legend([h1,h2,h3,h4,h5,h6],{'nomag1','nomag2','nomag3','mag1','mag2','mag3'});
xlim([0,27]);
% ylim([0,40]);

xvec = linspace(0.25,23.7,1000);
% xvec = linspace(0.25,26.4,1000);

Fnomag_mea = [interp1(data_nomag1_mea.depth_insertion(2:end-1),data_nomag1_mea.Fmag_smooth(2:end-1),xvec);...
              interp1(data_nomag2_mea.depth_insertion(2:end-1),data_nomag2_mea.Fmag_smooth(2:end-1),xvec);...
              interp1(data_nomag3_mea.depth_insertion(2:end-1),data_nomag3_mea.Fmag_smooth(2:end-1),xvec)];

Fmag = [interp1(data_mag1.depth_insertion(2:end-1),data_mag1.Fmag_smooth(2:end-1),xvec);...
        interp1(data_mag2.depth_insertion(2:end-1),data_mag2.Fmag_smooth(2:end-1),xvec);...
        interp1(data_mag3.depth_insertion(2:end-1),data_mag3.Fmag_smooth(2:end-1),xvec)];

Favg_nomag_mea = nanmean(Fnomag_mea,1);
std_nomag_mea = std(Fnomag_mea);

Favg_mag = nanmean(Fmag,1);
std_mag = std(Fmag);

colorsMat2 = distinguishable_colors(6);

% Plot the averages with no shifting
figure(2); grid on; hold on;
h1 = plot(xvec, Favg_nomag_mea, 'Color', colorsMat2(1,:), 'LineWidth',1,'LineStyle','--');
fill([xvec fliplr(xvec)],[Favg_nomag_mea + std_nomag_mea, fliplr(Favg_nomag_mea-std_nomag_mea)],...
    'b','FaceAlpha',0.2,'LineStyle','none');
h2 = plot(xvec, Favg_mag, 'Color', colorsMat2(2,:), 'LineWidth',1,'LineStyle','--');
fill([xvec fliplr(xvec)],[Favg_mag + std_mag, fliplr(Favg_mag-std_mag)],...
    'r','FaceAlpha',0.2,'LineStyle','none');

h4 = plot(xvec,-(Favg_nomag_mea-Favg_mag),'Color',colorsMat2(4,:),'LineWidth',4,'LineStyle',':');

legend([h1,h2,h4],{'Favg_{nomag_{mea}} N=3','Favg_{mag} N=3','Favgdiff N=3'});
xlabel('Insertion Depth [mm]'); ylabel('Force [mN]');