%% Data Analysis for Cadaver Trials
% Trevor Bruns and Katy Riojas
% Last Updated: 9/5/19

clear all; clc;

%% Inputs
addpath('classDef','functions','data');
load('avg_cal_slope.mat');
load('data\cadaver\cadaver_T_st_fixture.mat');
T_fixture_st = inv(T_st_fixture);
smooth_span = 0.06; % specify smoothing span

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

% Create data objects
data_nomag1_mea = MagneticGuidanceData(filepaths_nomag1_mea,avg_cal_slope);
data_nomag2_mea = MagneticGuidanceData(filepaths_nomag2_mea,avg_cal_slope);
data_nomag3_mea = MagneticGuidanceData(filepaths_nomag3_mea,avg_cal_slope);

data_mag1 = MagneticGuidanceData(filepaths_mag1,avg_cal_slope);
data_mag2 = MagneticGuidanceData(filepaths_mag2,avg_cal_slope);
data_mag3 = MagneticGuidanceData(filepaths_mag3,avg_cal_slope);

% Smooth Data
data_nomag1_mea = data_nomag1_mea.setSmoothSpan(smooth_span);
data_nomag2_mea = data_nomag2_mea.setSmoothSpan(smooth_span);
data_nomag3_mea = data_nomag3_mea.setSmoothSpan(smooth_span);
data_mag1 = data_mag1.setSmoothSpan(smooth_span);
data_mag2 = data_mag2.setSmoothSpan(smooth_span);
data_mag3 = data_mag3.setSmoothSpan(smooth_span);

% Rotate force data to align with ST frame
data_nomag1_mea = data_nomag1_mea.setT(T_fixture_st);
data_nomag2_mea = data_nomag2_mea.setT(T_fixture_st);
data_nomag3_mea = data_nomag3_mea.setT(T_fixture_st);
data_mag1 = data_mag1.setT(T_fixture_st);
data_mag2 = data_mag2.setT(T_fixture_st);
data_mag3 = data_mag3.setT(T_fixture_st);

%% Plot 1: XYZ raw and smoothed force data

alpha = 1; % line opacity
colorsMat = distinguishable_colors(6); % pull in a set of distinguishable colors
xyzColor = distinguishable_colors(3); % xyz colors

figure(1); clf(1);
sgtitle('Cadaver Insertion XYZ Forces');

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

%% Plot 2: Force Magnitude Comparison for Each Trial
figure(2); clf(2); hold on; grid on;
title('Guided vs Unguided Insertion')
xlabel('Insertion Depth [mm]')
ylabel('||Force|| [mN]')
xlim([0,27]);
% ylim([0,40]);

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
legend([h1,h2,h3,h4,h5,h6],{'nomag1','nomag2','nomag3','mag1','mag2','mag3'});

%% Plot 3: Averages
xvec = linspace(0.25,23.7,1000);

Fnomag_mea = [interp1(data_nomag1_mea.depth_insertion(2:end-1),...
                        data_nomag1_mea.Fmag_smooth(2:end-1),xvec);...
              interp1(data_nomag2_mea.depth_insertion(2:end-1),...
                        data_nomag2_mea.Fmag_smooth(2:end-1),xvec);...
              interp1(data_nomag3_mea.depth_insertion(2:end-1),...
                        data_nomag3_mea.Fmag_smooth(2:end-1),xvec)];

Fmag = [interp1(data_mag1.depth_insertion(2:end-1),...
                data_mag1.Fmag_smooth(2:end-1),xvec);...
        interp1(data_mag2.depth_insertion(2:end-1),...
                data_mag2.Fmag_smooth(2:end-1),xvec);...
        interp1(data_mag3.depth_insertion(2:end-1),...
                data_mag3.Fmag_smooth(2:end-1),xvec)];

% Compute Averages and Standard Deviations
Favg_nomag_mea = mean(Fnomag_mea,1);
std_nomag_mea = std(Fnomag_mea);

Favg_mag = mean(Fmag,1);
std_mag = std(Fmag);

% Plot the averages with no shifting
figure(3); clf(3); grid on; hold on;
xlabel('Insertion Depth [mm]'); ylabel('Average ||Force|| [mN]');
xlim([0,27]); ylim([0,120]);
xlim([0,xvec(end)]);
h1 = plot(xvec, Favg_nomag_mea, 'Color', 'r', 'LineWidth',1,'LineStyle','--');
fill([xvec fliplr(xvec)],[Favg_nomag_mea + std_nomag_mea, fliplr(Favg_nomag_mea-std_nomag_mea)],...
    'r','FaceAlpha',0.2,'LineStyle','none');
h2 = plot(xvec, Favg_mag, 'Color', 'b', 'LineWidth',1,'LineStyle','-');
fill([xvec fliplr(xvec)],[Favg_mag + std_mag, fliplr(Favg_mag-std_mag)],...
    'b','FaceAlpha',0.2,'LineStyle','none');
% h3 = plot(xvec,-(Favg_nomag_mea-Favg_mag),'Color','k','LineWidth',4,'LineStyle',':');

legend([h1,h2],{'No Magnetic Guidance N=3',...
    'Magnetic Guidance N=3'});%'Average Difference N=3'});

%% Plot 4: XYZ forces in ST frame
% Now plot in approximate st frame
figure(4); clf(4);
sgtitle('Cadaver XYZ Insertion Forces in ST Frame');

subplot(2,3,1); grid on; hold on;
ylabel('Unguided XYZ Forces [mN]');
xlim([0 28]); %ylim([-120,40]);
title('Trial 1');
% plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fx_smooth_st, 'Color', [xyzColor(1,:), 0.3*alpha]);
h1 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
h2 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
h3 = plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);

subplot(2,3,2); grid on; hold on;
xlim([0 28]); %ylim([-120,40]);
title('Trial 2');
% plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);

subplot(2,3,3); grid on; hold on;
xlim([0 28]); %ylim([-180,40]);
title('Trial 3');
% plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);

subplot(2,3,4); grid on; hold on;
xlim([0 28]); %ylim([-120,40]);
xlabel('Insertion Depth (mm)'); ylabel('Guided XYZ Forces [mN]');
% plot(data_mag1.depth_insertion, data_mag1.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_mag1.depth_insertion, data_mag1.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_mag1.depth_insertion, data_mag1.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_mag1.depth_insertion, data_mag1.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_mag1.depth_insertion, data_mag1.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_mag1.depth_insertion, data_mag1.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);

subplot(2,3,5); grid on; hold on;
xlim([0 28]); %ylim([-120,40]);
xlabel('Insertion Depth (mm)');
% plot(data_mag2.depth_insertion, data_mag2.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_mag2.depth_insertion, data_mag2.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_mag2.depth_insertion, data_mag2.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_mag2.depth_insertion, data_mag2.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_mag2.depth_insertion, data_mag2.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_mag2.depth_insertion, data_mag2.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);

subplot(2,3,6); grid on; hold on;
xlim([0 28]); %ylim([-120,40]);
xlabel('Insertion Depth (mm)');
% plot(data_mag3.depth_insertion, data_mag3.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_mag3.depth_insertion, data_mag3.Fx_smooth_st, 'Color', [xyzColor(1,:), alpha]);
% plot(data_mag3.depth_insertion, data_mag3.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_mag3.depth_insertion, data_mag3.Fy_smooth_st, 'Color', [xyzColor(2,:), alpha]);
% plot(data_mag3.depth_insertion, data_mag3.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
plot(data_mag3.depth_insertion, data_mag3.Fz_smooth_st, 'Color', [xyzColor(3,:), alpha]);
legend([h1,h2,h3],'-Fz','Fx','-Fy');

%% Save Raw Force Magnitudes
c_ug_Fmag = [data_nomag1_mea.Fmag;data_nomag2_mea.Fmag;data_nomag3_mea.Fmag];
c_g_Fmag = [data_mag1.Fmag;data_mag2.Fmag;data_mag3.Fmag];

save('data\cadaver\c_ug_Fmag.mat','c_ug_Fmag');
save('data\cadaver\c_g_Fmag.mat','c_g_Fmag');