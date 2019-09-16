%Data Analysis for Cadaver Trials
% Trevor Bruns and Katy Riojas
% Last Updated: 9/5/19

%% Inputs
clear all; clc;
addpath('classDef','functions','data');
load('avg_cal_slopes.mat');
load('data\cadaver\cadaver_T_st_fixture.mat');
T_fixture_st = inv(T_st_fixture);

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
data_nomag1_mea = MagneticGuidanceData(filepaths_nomag1_mea, cal_slopes);
data_nomag2_mea = MagneticGuidanceData(filepaths_nomag2_mea, cal_slopes);
data_nomag3_mea = MagneticGuidanceData(filepaths_nomag3_mea, cal_slopes);

data_mag1 = MagneticGuidanceData(filepaths_mag1, cal_slopes);
data_mag2 = MagneticGuidanceData(filepaths_mag2, cal_slopes);
data_mag3 = MagneticGuidanceData(filepaths_mag3, cal_slopes);

% Set smoothing span (proportion of data points)
force_smooth_span = 40; % [# samples]

% Smooth Data
data_nomag1_mea.smooth_span = force_smooth_span;
data_nomag2_mea.smooth_span = force_smooth_span;
data_nomag3_mea.smooth_span = force_smooth_span;

data_mag1.smooth_span = force_smooth_span;
data_mag2.smooth_span = force_smooth_span;
data_mag3.smooth_span = force_smooth_span;

%% Trim to 125mN Max
Fthresh = 125;
nomag1_endi = find(data_nomag1_mea.Fmag>Fthresh,1);
if isempty(nomag1_endi), nomag1_endi = length(data_nomag1_mea.Fmag); end
nomag2_endi = find(data_nomag2_mea.Fmag>Fthresh,1);
if isempty(nomag2_endi), nomag2_endi = length(data_nomag2_mea.Fmag); end
nomag3_endi = find(data_nomag3_mea.Fmag>Fthresh,1);
if isempty(nomag3_endi), nomag3_endi = length(data_nomag3_mea.Fmag); end
mag1_endi = find(data_mag1.Fmag>Fthresh,1);
if isempty(mag1_endi), mag1_endi = length(data_mag1.Fmag); end
mag2_endi = find(data_mag2.Fmag>Fthresh,1);
if isempty(mag2_endi), mag2_endi = length(data_mag2.Fmag); end
mag3_endi = find(data_mag3.Fmag>Fthresh,1);
if isempty(mag3_endi), mag3_endi = length(data_mag3.Fmag); end

%% Pull in Plot variables
alpha = 1; % line opacity
colorsMat = distinguishable_colors(6); % pull in a set of distinguishable colors
xyzColor = distinguishable_colors(3); % xyz colors
maxY = max([data_nomag1_mea.Fmag;data_nomag2_mea.Fmag;...
    data_nomag3_mea.Fmag;data_mag1.Fmag;data_mag2.Fmag;data_mag3.Fmag]);
plotXYZdata = 0; %toggle whether to do plot3

%% Plot 2: Force Magnitude Comparison for Each Trial
figure(1); clf(1);
subplot(2,3,1); grid on; hold on;
sgtitle('Raw and Smoothed Fmag Cadaver Trials');
ylabel('Unguided ||Force|| [mN]');
ylim([0 maxY]);
% xlabel('Insertion Depth [mm]'); 
title('Trial 1');
plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fmag, 'Color', [colorsMat(1,:),  0.3*alpha], 'LineWidth',1);
plot(data_nomag1_mea.depth_insertion, data_nomag1_mea.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth',1);
scatter(data_nomag1_mea.depth_insertion(nomag1_endi),data_nomag1_mea.Fmag(nomag1_endi),10,'filled','r');

subplot(2,3,2); grid on; hold on;
% xlabel('Insertion Depth [mm]'); 
title('Trial 2');
ylim([0 maxY]);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fmag, 'Color', [colorsMat(1,:), 0.3*alpha], 'LineWidth',1);
plot(data_nomag2_mea.depth_insertion, data_nomag2_mea.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth',1);
scatter(data_nomag2_mea.depth_insertion(nomag2_endi),data_nomag2_mea.Fmag(nomag2_endi),10,'filled','r');

subplot(2,3,3); grid on; hold on;
% xlabel('Insertion Depth [mm]'); 
title('Trial 3');
ylim([0 maxY]);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fmag, 'Color', [colorsMat(1,:),  0.3*alpha], 'LineWidth',1);
plot(data_nomag3_mea.depth_insertion, data_nomag3_mea.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth',1);
scatter(data_nomag3_mea.depth_insertion(nomag3_endi),data_nomag3_mea.Fmag(nomag3_endi),10,'filled','r');

subplot(2,3,4); grid on; hold on;
xlabel('Insertion Depth [mm]'); title('Trial 1');
ylabel('Guided ||Force|| [mN]');
ylim([0 maxY]);
plot(data_mag1.depth_insertion, data_mag1.Fmag, 'Color', [colorsMat(1,:), 0.3*alpha], 'LineWidth', 1);
h4 = plot(data_mag1.depth_insertion, data_mag1.Fmag_smooth, 'Color',colorsMat(1,:), 'LineWidth', 2);
scatter(data_mag1.depth_insertion(mag1_endi),data_mag1.Fmag(mag1_endi),10,'filled','r');

subplot(2,3,5); grid on; hold on;
xlabel('Insertion Depth [mm]'); title('Trial 2');
ylim([0 maxY]);
plot(data_mag2.depth_insertion, data_mag2.Fmag, 'Color', [colorsMat(1,:), 0.3*alpha], 'LineWidth', 1);
h5 = plot(data_mag2.depth_insertion, data_mag2.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth', 2);
scatter(data_mag2.depth_insertion(mag2_endi),data_mag2.Fmag(mag2_endi),10,'filled','r');

subplot(2,3,6); grid on; hold on;
xlabel('Insertion Depth [mm]'); title('Trial 3');
ylim([0 maxY]);
plot(data_mag3.depth_insertion, data_mag3.Fmag, 'Color', [colorsMat(1,:), 0.3*alpha], 'LineWidth', 1);
h6 = plot(data_mag3.depth_insertion, data_mag3.Fmag_smooth, 'Color', colorsMat(1,:), 'LineWidth', 2);
scatter(data_mag3.depth_insertion(mag3_endi),data_mag3.Fmag(mag3_endi),10,'filled','r');

%% Plot 2: Averages
maxX = min([data_nomag1_mea.depth_insertion(nomag1_endi);...
            data_nomag2_mea.depth_insertion(nomag2_endi);...
            data_nomag3_mea.depth_insertion(nomag3_endi);...
            data_mag1.depth_insertion(mag1_endi);...
            data_mag2.depth_insertion(mag2_endi);...
            data_mag3.depth_insertion(mag3_endi)]);
            
maxY = min([data_nomag1_mea.Fmag_smooth(nomag1_endi);...
            data_nomag2_mea.Fmag_smooth(nomag2_endi);...
            data_nomag3_mea.Fmag_smooth(nomag3_endi);...
            data_mag1.Fmag_smooth(mag1_endi);...
            data_mag2.Fmag_smooth(mag2_endi);...
            data_mag3.Fmag_smooth(mag3_endi)]);

xvec = linspace(0.02,maxX,1000);

Fnomag_mea = [interp1(data_nomag1_mea.depth_insertion(1:nomag1_endi),...
                        data_nomag1_mea.Fmag_smooth(1:nomag1_endi),xvec);...
              interp1(data_nomag2_mea.depth_insertion(1:nomag2_endi),...
                        data_nomag2_mea.Fmag_smooth(1:nomag2_endi),xvec);...
              interp1(data_nomag3_mea.depth_insertion(1:nomag3_endi),...
                        data_nomag3_mea.Fmag_smooth(1:nomag3_endi),xvec)];

Fmag = [interp1(data_mag1.depth_insertion(1:mag1_endi),...
                data_mag1.Fmag_smooth(1:mag1_endi),xvec);...
        interp1(data_mag2.depth_insertion(1:mag2_endi),...
                data_mag2.Fmag_smooth(1:mag2_endi),xvec);...
        interp1(data_mag3.depth_insertion(1:mag3_endi),...
                data_mag3.Fmag_smooth(1:mag3_endi),xvec)];

% Compute Averages and Standard Deviations
Favg_nomag_mea = mean(Fnomag_mea,1);
std_nomag_mea = std(Fnomag_mea);

Favg_mag = mean(Fmag,1);
std_mag = std(Fmag);

% Plot the averages with no shifting
figure(2); clf(2); grid on; hold on;
xlabel('Insertion Depth [mm]'); ylabel('Average ||Force|| [mN]');
xlim([0,maxX]); %ylim([0,maxY]);
ylim([0 120]);
h1 = plot(xvec, Favg_nomag_mea, 'Color', 'r', 'LineWidth',1,'LineStyle','--');
fill([xvec fliplr(xvec)],[Favg_nomag_mea + std_nomag_mea, fliplr(Favg_nomag_mea-std_nomag_mea)],...
    'r','FaceAlpha',0.2,'LineStyle','none');
h2 = plot(xvec, Favg_mag, 'Color', 'b', 'LineWidth',1,'LineStyle','-');
fill([xvec fliplr(xvec)],[Favg_mag + std_mag, fliplr(Favg_mag-std_mag)],...
    'b','FaceAlpha',0.2,'LineStyle','none');

legend([h1,h2],{'Case 2','Case 3'});

%% Save Raw Force Magnitudes
c_ug_Fmag = [data_nomag1_mea.Fmag(1:nomag1_endi);...
             data_nomag2_mea.Fmag(1:nomag2_endi);...
             data_nomag3_mea.Fmag(1:nomag3_endi)];
c_g_Fmag = [data_mag1.Fmag(1:mag1_endi);...
            data_mag2.Fmag(1:mag2_endi);...
            data_mag3.Fmag(1:mag3_endi)];

save('data\cadaver\c_ug_Fmag.mat','c_ug_Fmag');
save('data\cadaver\c_g_Fmag.mat','c_g_Fmag');

%% Plot 3: XYZ raw and smoothed force data
if plotXYZdata
    figure(3); clf(3);
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
end