%% import force and smaract data
clear all; close all; clc;

addpath('classDef','functions','data');
load('avg_cal_slope.mat');

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

% create data objects
data_nomag1_mea = MagneticGuidanceData(filepaths_nomag1_mea);
data_nomag2_mea = MagneticGuidanceData(filepaths_nomag2_mea);
data_nomag3_mea = MagneticGuidanceData(filepaths_nomag3_mea);
data_nomag4_mea = MagneticGuidanceData(filepaths_nomag4_mea);

data_mag1 = MagneticGuidanceData(filepaths_mag1);
data_mag2 = MagneticGuidanceData(filepaths_mag2);
data_mag3 = MagneticGuidanceData(filepaths_mag3);
data_mag4 = MagneticGuidanceData(filepaths_mag4);

data_nomag1_ea = MagneticGuidanceData(filepaths_nomag1_ea);
data_nomag2_ea = MagneticGuidanceData(filepaths_nomag2_ea);
data_nomag3_ea = MagneticGuidanceData(filepaths_nomag3_ea);
data_nomag4_ea = MagneticGuidanceData(filepaths_nomag4_ea);

data_nomag1_mea = data_nomag1_mea.setCalSlope(avg_cal_slope);
data_nomag2_mea = data_nomag2_mea.setCalSlope(avg_cal_slope);
data_nomag3_mea = data_nomag3_mea.setCalSlope(avg_cal_slope);
data_nomag4_mea = data_nomag4_mea.setCalSlope(avg_cal_slope);

data_mag1 = data_mag1.setCalSlope(avg_cal_slope);
data_mag2 = data_mag2.setCalSlope(avg_cal_slope);
data_mag3 = data_mag3.setCalSlope(avg_cal_slope);
data_mag4 = data_mag4.setCalSlope(avg_cal_slope);

data_nomag1_ea = data_nomag1_ea.setCalSlope(avg_cal_slope);
data_nomag2_ea = data_nomag2_ea.setCalSlope(avg_cal_slope);
data_nomag3_ea = data_nomag3_ea.setCalSlope(avg_cal_slope);
data_nomag4_ea = data_nomag4_ea.setCalSlope(avg_cal_slope);

% set smoothing span (proportion of data points)
smooth_span = 0.06;

data_nomag1_mea = data_nomag1_mea.setSmoothSpan(smooth_span);
data_nomag2_mea = data_nomag2_mea.setSmoothSpan(smooth_span);
data_nomag3_mea = data_nomag3_mea.setSmoothSpan(smooth_span);
data_nomag4_mea = data_nomag4_mea.setSmoothSpan(smooth_span);

data_mag1 = data_mag1.setSmoothSpan(smooth_span);
data_mag2 = data_mag2.setSmoothSpan(smooth_span);
data_mag3 = data_mag3.setSmoothSpan(smooth_span);
data_mag4 = data_mag4.setSmoothSpan(smooth_span);

data_nomag1_ea = data_nomag1_ea.setSmoothSpan(smooth_span);
data_nomag2_ea = data_nomag2_ea.setSmoothSpan(smooth_span);
data_nomag3_ea = data_nomag3_ea.setSmoothSpan(smooth_span);
data_nomag4_ea = data_nomag4_ea.setSmoothSpan(smooth_span);

%% plot
PlotColors; % load custom colors
alpha = 1; % reduce transparency of unguided plot lines

colorsMat = distinguishable_colors(12);

% Figure: force magnitude
figure(2); hold on; grid on;
plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fmag_cal, 'Color', [colorsMat(1,:),  0.3*alpha], 'LineWidth',1);
h1 = plot(data_nomag1_ea.depth_insertion, data_nomag1_ea.Fmag_smooth_cal, 'Color', colorsMat(1,:), 'LineWidth',3,'LineStyle',':');
% 
plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fmag_cal, 'Color', colorsMat(2,:), 'LineWidth',1);
h2 = plot(data_nomag2_ea.depth_insertion, data_nomag2_ea.Fmag_smooth_cal, 'Color', colorsMat(2,:), 'LineWidth',3,'LineStyle',':');

plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fmag_cal, 'Color', [colorsMat(3,:),  0.3*alpha], 'LineWidth',1);
h3 = plot(data_nomag3_ea.depth_insertion, data_nomag3_ea.Fmag_smooth_cal, 'Color', colorsMat(3,:), 'LineWidth',3,'LineStyle',':');

plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fmag_cal, 'Color', [colorsMat(4,:),  0.3*alpha], 'LineWidth',1);
h4 = plot(data_nomag4_ea.depth_insertion, data_nomag4_ea.Fmag_smooth_cal, 'Color', colorsMat(4,:), 'LineWidth',3,'LineStyle',':');

plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fmag_cal, 'Color', [colorsMat(5,:),  0.3*alpha], 'LineWidth',1);
h5 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fmag_smooth_cal, 'Color', colorsMat(5,:), 'LineWidth',1,'LineStyle','--');
% 
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fmag_cal, 'Color', colorsMat(6,:), 'LineWidth',1);
h6 = plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fmag_smooth_cal, 'Color', colorsMat(6,:), 'LineWidth',1,'LineStyle','--');

plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fmag_cal, 'Color', [colorsMat(7,:),  0.3*alpha], 'LineWidth',1);
h7 = plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fmag_smooth_cal, 'Color', colorsMat(7,:), 'LineWidth',1,'LineStyle','--');

plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fmag_cal, 'Color', [colorsMat(8,:),  0.3*alpha], 'LineWidth',1);
h8 = plot(data_nomag4_mea.depth_insertion, data_nomag4_mea.Fmag_smooth_cal, 'Color', colorsMat(8,:), 'LineWidth',1,'LineStyle','--');

plot(data_mag1.depth_insertion, data_mag1.Fmag_cal, 'Color', [colorsMat(9,:), 0.3*alpha], 'LineWidth', 1);
h9 = plot(data_mag1.depth_insertion, data_mag1.Fmag_smooth_cal, 'Color',colorsMat(9,:), 'LineWidth', 2);

plot(data_mag2.depth_insertion, data_mag2.Fmag_cal, 'Color', [colorsMat(10,:), 0.3*alpha], 'LineWidth', 1);
h10 = plot(data_mag2.depth_insertion, data_mag2.Fmag_smooth_cal, 'Color', colorsMat(10,:), 'LineWidth', 2);

plot(data_mag3.depth_insertion, data_mag3.Fmag_cal, 'Color', [colorsMat(11,:), 0.3*alpha], 'LineWidth', 1);
h11 = plot(data_mag3.depth_insertion, data_mag3.Fmag_smooth_cal, 'Color', colorsMat(11,:), 'LineWidth', 2);

plot(data_mag4.depth_insertion, data_mag4.Fmag_cal, 'Color', [colorsMat(12,:), 0.3*alpha], 'LineWidth', 1);
h12 = plot(data_mag4.depth_insertion, data_mag4.Fmag_smooth_cal, 'Color', colorsMat(12,:), 'LineWidth', 2);

title('Guided vs Unguided Insertion')
xlabel('Insertion Depth [mm]')
ylabel('||Force|| [mN]')
% legend([h1,h4],{'nomag1','mag1'});
legend([h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12],{'nomag1-ea',...
    'nomag2-ea','nomag3-ea','nomag4-ea','nomag1-mea',...
    'nomag2-mea','nomag3-mea','nomag4-mea','mag1','mag2','mag3','mag4'});


xvec = linspace(0.25,23.7,1000);
% xvec = linspace(0.25,26.4,1000);

Fnomag_mea = [interp1(data_nomag1_mea.depth_insertion(2:end-1),data_nomag1_mea.Fmag_smooth_cal(2:end-1),xvec);...
              interp1(data_nomag2_mea.depth_insertion(2:end-1),data_nomag2_mea.Fmag_smooth_cal(2:end-1),xvec);...
              interp1(data_nomag3_mea.depth_insertion(2:end-1),data_nomag3_mea.Fmag_smooth_cal(2:end-1),xvec);...
              interp1(data_nomag4_mea.depth_insertion(2:end-1),data_nomag4_mea.Fmag_smooth_cal(2:end-1),xvec)];
% 
Fmag = [interp1(data_mag1.depth_insertion(2:end-1),data_mag1.Fmag_smooth_cal(2:end-1),xvec);...
        interp1(data_mag2.depth_insertion(2:end-1),data_mag2.Fmag_smooth_cal(2:end-1),xvec);...
        interp1(data_mag3.depth_insertion(2:end-1),data_mag3.Fmag_smooth_cal(2:end-1),xvec);...
        interp1(data_mag4.depth_insertion(2:end-1),data_mag4.Fmag_smooth_cal(2:end-1),xvec)];
    
Fnomag_ea = [interp1(data_nomag1_ea.depth_insertion(2:end-1),data_nomag1_ea.Fmag_smooth_cal(2:end-1),xvec);...
             interp1(data_nomag2_ea.depth_insertion(2:end-1),data_nomag2_ea.Fmag_smooth_cal(2:end-1),xvec);...
             interp1(data_nomag3_ea.depth_insertion(2:end-1),data_nomag3_ea.Fmag_smooth_cal(2:end-1),xvec);...
             interp1(data_nomag4_ea.depth_insertion(2:end-1),data_nomag4_ea.Fmag_smooth_cal(2:end-1),xvec)];

Favg_nomag_mea = mean(Fnomag_mea,1);
std_nomag_mea = std(Fnomag_mea);

Favg_mag = mean(Fmag,1);
std_mag = std(Fmag);

Favg_nomag_ea = nanmean(Fnomag_ea,1);
std_nomag_ea = std(Fnomag_ea);

colorsMat2 = distinguishable_colors(6);

% Plot the averages with no shifting
figure(3); grid on; hold on;
xlim([0,27]);
ylim([-40,120]);
h1 = plot(xvec, Favg_nomag_mea, 'Color', colorsMat2(1,:), 'LineWidth',1,'LineStyle','--');
fill([xvec fliplr(xvec)],[Favg_nomag_mea + std_nomag_mea, fliplr(Favg_nomag_mea-std_nomag_mea)],...
    'b','FaceAlpha',0.2,'LineStyle','none');

h2 = plot(xvec, Favg_mag, 'Color', colorsMat2(2,:), 'LineWidth',2,'LineStyle','-');
fill([xvec fliplr(xvec)],[Favg_mag + std_mag, fliplr(Favg_mag-std_mag)],...
    'r','FaceAlpha',0.2,'LineStyle','none');

% h3 = plot(xvec, Favg_nomag_ea, 'Color', colorsMat2(3,:), 'LineWidth',3,'LineStyle',':');
% fill([xvec fliplr(xvec)],[Favg_nomag_ea + std_nomag_ea, fliplr(Favg_nomag_ea-std_nomag_ea)],...
%     'g','FaceAlpha',0.2,'LineStyle','none');

h4 = plot(xvec,-(Favg_nomag_mea-Favg_mag),'Color',colorsMat2(4,:),'LineWidth',4,'LineStyle',':');

% Plot differences
% h5 = plot(xvec,(Favg_nomag_mea-Favg_nomag_ea),'Color',colorsMat2(5,:),'LineWidth',4,'LineStyle',':');
% h6 = plot(xvec,-(Favg_nomag_ea-Favg_mag),'Color',colorsMat2(6,:),'LineWidth',4,'LineStyle',':');
legend([h1,h2,h4],'No Magnetic Guidance N=4','Magnetic Guidance N=4','Average Difference N=4');
% legend([h1,h2,h3,h4,h5,h6],{'Favg_{nomag_mea}','Favg_{mag}',...
%     'Favg_{nomag_ea}','Favgdiff_{nomagMEA2mag}','Favgdiff_{nomagMEA2nomagEA}',...
%     'Favgdiff_{nomagEA2mag}'});
xlabel('Insertion Depth [mm]'); ylabel('Average ||Force|| [mN]');