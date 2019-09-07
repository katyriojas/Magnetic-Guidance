%% import force and smaract data
% clear all; close all; clc;

addpath('classDef','functions','data');
load('avg_cal_slope.mat');
% Specify filepaths for .csv files of data
base_path = 'data\cadaver\';

filepaths_nomag1_mea.smaract = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial1_smaract.csv');
filepaths_nomag1_mea.force   = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial1_force.csv');

filepaths_nomag2_mea.smaract = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial2_smaract.csv');
filepaths_nomag2_mea.force   = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial2_force.csv');

filepaths_nomag3_mea.smaract = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial3_smaract.csv');
filepaths_nomag3_mea.force   = fullfile(base_path, 'UG-MEA\cadaver_ug_mea2_trial3_force.csv');

filepaths_mag1.smaract   = fullfile(base_path, 'guided\cadaver_g_mea2_trial1_smaract.csv');
filepaths_mag1.force     = fullfile(base_path, 'guided\cadaver_g_mea2_trial1_force.csv');

filepaths_mag2.smaract = fullfile(base_path, 'guided\cadaver_g_mea2_trial2_smaract.csv');
filepaths_mag2.force   = fullfile(base_path, 'guided\cadaver_g_mea2_trial2_force.csv');

filepaths_mag3.smaract = fullfile(base_path, 'guided\cadaver_g_mea2_trial3_smaract.csv');
filepaths_mag3.force   = fullfile(base_path, 'guided\cadaver_g_mea2_trial3_force.csv');

%% Create data objects
data_nomag1_mea = MagneticGuidanceData(filepaths_nomag1_mea, avg_cal_slope);
% data_nomag1_mea = MagneticGuidanceData(filepaths_nomag1_mea);
data_nomag2_mea = MagneticGuidanceData(filepaths_nomag2_mea, avg_cal_slope);
data_nomag3_mea = MagneticGuidanceData(filepaths_nomag3_mea, avg_cal_slope);

data_mag1 = MagneticGuidanceData(filepaths_mag1, avg_cal_slope);
data_mag2 = MagneticGuidanceData(filepaths_mag2, avg_cal_slope);
data_mag3 = MagneticGuidanceData(filepaths_mag3, avg_cal_slope);


% Set smoothing span (proportion of data points)
smooth_span = 0.1;

data_nomag1_mea = data_nomag1_mea.setSmoothSpan(smooth_span);
% data_nomag1_mea.setSmoothSpan(smooth_span);
data_nomag2_mea = data_nomag2_mea.setSmoothSpan(smooth_span);
data_nomag3_mea = data_nomag3_mea.setSmoothSpan(smooth_span);
data_mag1 = data_mag1.setSmoothSpan(smooth_span);
data_mag2 = data_mag2.setSmoothSpan(smooth_span);
data_mag3 = data_mag3.setSmoothSpan(smooth_span);



%% Plot
alpha = 1; % reduce transparency of unguided plot lines
colorsMat = distinguishable_colors(12); % pull in a set of distinguishable colors
xyzColor = distinguishable_colors(3); % xyz colors

figure(1); clf(1);
sgtitle('Cadaver Insertions');

subplot(2,3,1); grid on; hold on;
ylabel('Unguided XYZ Forces [mN]');
xlim([0 28]); ylim([-120,40]);
title('Trial 1');
plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
h1 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
h2 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
h3 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(2,3,2); grid on; hold on;
xlim([0 28]); ylim([-120,40]);
title('Trial 2');
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(2,3,3); grid on; hold on;
xlim([0 28]); ylim([-120,40]);
title('Trial 3');
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(2,3,4); grid on; hold on;
xlim([0 28]); ylim([-120,40]);
xlabel('Insertion Depth (mm)'); ylabel('Guided XYZ Forces [mN]');
plot(data_mag1.depth_insertion, data_mag1.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_mag1.depth_insertion, data_mag1.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
plot(data_mag1.depth_insertion, data_mag1.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_mag1.depth_insertion, data_mag1.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
plot(data_mag1.depth_insertion, data_mag1.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_mag1.depth_insertion, data_mag1.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(2,3,5); grid on; hold on;
xlim([0 28]); ylim([-120,40]);
xlabel('Insertion Depth (mm)');
plot(data_mag2.depth_insertion, data_mag2.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_mag2.depth_insertion, data_mag2.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
plot(data_mag2.depth_insertion, data_mag2.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_mag2.depth_insertion, data_mag2.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
plot(data_mag2.depth_insertion, data_mag2.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_mag2.depth_insertion, data_mag2.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
legend([h1,h2,h3],'Fx','Fy','Fz');

subplot(2,3,6); grid on; hold on;
xlim([0 28]); ylim([-120,40]);
xlabel('Insertion Depth (mm)');
plot(data_mag3.depth_insertion, data_mag3.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_mag3.depth_insertion, data_mag3.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
plot(data_mag3.depth_insertion, data_mag3.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_mag3.depth_insertion, data_mag3.Fy_smooth, 'Color', [xyzColor(2,:), alpha]);
plot(data_mag3.depth_insertion, data_mag3.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_mag3.depth_insertion, data_mag3.Fz_smooth, 'Color', [xyzColor(3,:), alpha]);
legend([h1,h2,h3],'Fx','Fy','Fz');

% Create force magnitude figure
figure(2); clf(2); hold on; grid on;

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

Favg_nomag_mea = mean(Fnomag_mea,1);
std_nomag_mea = std(Fnomag_mea);

Favg_mag = mean(Fmag,1);
std_mag = std(Fmag);

colorsMat2 = distinguishable_colors(6);

% Plot the averages with no shifting
figure(3); clf(3); grid on; hold on;
xlim([0,27]);
ylim([-40,120]);
h1 = plot(xvec, Favg_nomag_mea, 'Color', colorsMat2(1,:), 'LineWidth',1,'LineStyle','--');
fill([xvec fliplr(xvec)],[Favg_nomag_mea + std_nomag_mea, fliplr(Favg_nomag_mea-std_nomag_mea)],...
    'b','FaceAlpha',0.2,'LineStyle','none');
h2 = plot(xvec, Favg_mag, 'Color', colorsMat2(2,:), 'LineWidth',1,'LineStyle','--');
fill([xvec fliplr(xvec)],[Favg_mag + std_mag, fliplr(Favg_mag-std_mag)],...
    'r','FaceAlpha',0.2,'LineStyle','none');

h4 = plot(xvec,-(Favg_nomag_mea-Favg_mag),'Color',colorsMat2(4,:),'LineWidth',4,'LineStyle',':');

legend([h1,h2,h4],{'No Magnetic Guidance N=3','Magnetic Guidance N=3','Average Difference N=3'});
xlabel('Insertion Depth [mm]'); ylabel('Average ||Force|| [mN]');