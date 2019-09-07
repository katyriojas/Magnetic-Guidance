%% Import Force and Smaract data
% clear all; close all; clc;

addpath('functions','classDef','data');

% load force sensor calibration slopes
load('avg_cal_slope.mat');

% Filepaths to CSVs that were exported from ROS bags
base_path = 'data';

filepaths_p_manual.force1 = fullfile(base_path, 'phantom\manual\phantom_manual_ea_trial1_force.csv');
filepaths_p_manual.force2 = fullfile(base_path, 'phantom\manual\phantom_manual_ea_trial2_force.csv');
filepaths_p_manual.force3 = fullfile(base_path, 'phantom\manual\phantom_manual_ea_trial3_force.csv');
filepaths_p_manual.force4 = fullfile(base_path, 'phantom\manual\phantom_manual_ea2_trial4_force.csv');

filepaths_c_manual.force1 = fullfile(base_path, 'cadaver\manual\cadaver_manual_ea3_trial1_force.csv');
filepaths_c_manual.force2 = fullfile(base_path, 'cadaver\manual\cadaver_manual_ea3_trial2_force.csv');
filepaths_c_manual.force3 = fullfile(base_path, 'cadaver\manual\cadaver_manual_ea3_trial3_force.csv');

%% Create data objects
data_pman1 = Nano17Data(filepaths_p_manual.force1, cal_slopes);
data_pman2 = Nano17Data(filepaths_p_manual.force2, cal_slopes);
data_pman3 = Nano17Data(filepaths_p_manual.force3, cal_slopes);
data_pman4 = Nano17Data(filepaths_p_manual.force4, cal_slopes);

data_cman1 = Nano17Data(filepaths_c_manual.force1, cal_slopes);
data_cman2 = Nano17Data(filepaths_c_manual.force2, cal_slopes);
data_cman3 = Nano17Data(filepaths_c_manual.force3, cal_slopes);

%% Set smoothing
smooth_span = 40; % [# samples]
data_pman1.smooth_span = smooth_span;
data_pman2.smooth_span = smooth_span;
data_pman3.smooth_span = smooth_span;
data_pman4.smooth_span = smooth_span;

data_cman1.smooth_span = smooth_span;
data_cman2.smooth_span = smooth_span;
data_cman3.smooth_span = smooth_span;

%% Plot the manual data
alpha = 1; % reduce transparency of unguided plot lines
xyzColor = distinguishable_colors(3);
cMat = distinguishable_colors(7);

figure(1); 
subplot(1,4,1); grid on; hold on;
xlabel('Time (s)'); ylabel('Force [mN]');
sgtitle('Manual Insertions in Phantoms');

% manual, no smoothing
plot(data_pman1.time,data_pman1.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
% plot(data_pman1.time-data_pman1.time(1),data_pman1.Fx_smooth, 'Color', [xyzColor(1,:), alpha]);
plot(data_pman1.time,data_pman1.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_pman1.time,data_pman1.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);

subplot(1,4,2); grid on; hold on;
xlabel('Time (s)'); 
plot(data_pman2.time,data_pman2.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_pman2.time,data_pman2.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_pman2.time,data_pman2.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);

subplot(1,4,3); grid on; hold on;
xlabel('Time (s)'); 
plot(data_pman3.time,data_pman3.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_pman3.time,data_pman3.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_pman3.time,data_pman3.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);

subplot(1,4,4); grid on; hold on;
xlabel('Time (s)'); 
plot(data_pman4.time,data_pman4.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_pman4.time,data_pman4.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_pman4.time,data_pman4.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);
legend('Fx','Fy','Fz');

figure(2);
subplot(1,3,1); grid on; hold on;
xlabel('Time (s)'); ylabel('Force [mN]');
sgtitle('Manual Insertions in a Cadaver');

plot(data_cman1.time,data_cman1.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_cman1.time,data_cman1.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_cman1.time,data_cman1.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);

subplot(1,3,2); grid on; hold on;
plot(data_cman2.time,data_cman2.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_cman2.time,data_cman2.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_cman2.time,data_cman2.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);

subplot(1,3,3); grid on; hold on;
plot(data_cman3.time,data_cman3.Fx, 'Color', [xyzColor(1,:), 0.3*alpha]);
plot(data_cman3.time,data_cman3.Fy, 'Color', [xyzColor(2,:), 0.3*alpha]);
plot(data_cman3.time,data_cman3.Fz, 'Color', [xyzColor(3,:), 0.3*alpha]);

legend('Fx','Fy','Fz');

figure(3); 
subplot(1,2,1); grid on; hold on;
xlabel('Time (s)'); ylabel('Force [mN]');
title('Magnitude of Manual Insertions in Phantoms');

plot(data_pman1.time, data_pman1.Fmag, 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h1 = plot(data_pman1.time, data_pman1.Fmag_smooth, 'Color', [cMat(1,:),  alpha], 'LineWidth',1);

plot(data_pman2.time, data_pman2.Fmag, 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h2 = plot(data_pman2.time, data_pman2.Fmag_smooth, 'Color', [cMat(2,:),  0.3*alpha], 'LineWidth',1);

plot(data_pman3.time, data_pman3.Fmag, 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h3 = plot(data_pman3.time, data_pman3.Fmag_smooth, 'Color', [cMat(3,:),  0.3*alpha], 'LineWidth',1);

plot(data_pman4.time, data_pman4.Fmag, 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h4 = plot(data_pman4.time, data_pman4.Fmag_smooth, 'Color', [cMat(4,:),  0.3*alpha], 'LineWidth',1);

legend([h1,h2,h3,h4],{'Trial 1','Trial 2','Trial 3','Trial 4'});

subplot(1,2,2); grid on; hold on;
xlabel('Time (s)'); ylabel('Force [mN]');
title('Magnitude of Manual Insertions in Cadavers');

% plot(data_cman1.time, data_cman1.Fmag, 'Color', [cMat(5,:),  0.3*alpha], 'LineWidth',0.1);
h5 = plot(data_cman1.time, data_cman1.Fmag_smooth, 'Color', [cMat(5,:),  0.3*alpha], 'LineWidth',1);

% plot(data_cman2.time, data_cman2.Fmag, 'Color', [cMat(6,:),  0.3*alpha], 'LineWidth',0.1);
h6 = plot(data_cman2.time, data_cman2.Fmag_smooth, 'Color', [cMat(6,:),  0.3*alpha], 'LineWidth',1);

% plot(data_cman3.time, data_cman3.Fmag, 'Color', [cMat(7,:),  0.3*alpha], 'LineWidth',0.1);
h7= plot(data_cman3.time, data_cman3.Fmag_smooth, 'Color', [cMat(7,:),  0.3*alpha], 'LineWidth',1);

legend([h5,h6,h7],{'Trial 1','Trial 2','Trial 3'});

% xvec = linspace(0,100,1000);
% 
% Fmag_smooth_c = [interp1(10^-9*(data_cman1.time-data_cman1.time(1)),cman1_Fmag_smooth,xvec);...
%                  interp1(10^-9*(data_cman2.time-data_cman2.time(1)),cman2_Fmag_smooth,xvec);...
%                  interp1(10^-9*(data_cman3.time-data_cman3.time(1)),cman3_Fmag_smooth,xvec)];
% 
% Favg_manual_c = nanmean(Fmag_smooth_c,1);
% std_manual_c = nanstd(Fmag_smooth_c);
% 
% figure(4); 
% grid on; hold on; xlabel('Time (s)'); ylabel('Average ||Force||');
% title('Average Insertion Force in Cadaver vs. Time');
% plot(xvec,Favg_manual_c,'Color',cMat(1,:));
% fill([xvec fliplr(xvec)],[Favg_manual_c + std_manual_c, fliplr(Favg_manual_c-std_manual_c)],...
%     'b','FaceAlpha',0.2,'LineStyle','none');
