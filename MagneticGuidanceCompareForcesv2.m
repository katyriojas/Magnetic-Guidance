%% Import Force and Smaract data
clear all; close all; clc;

addpath('C:\Users\riojaske\Documents\magsteer\PreopPlanning\functions');
addpath('C:\Users\riojaske\Documents\magsteer\PreopPlanning\classDef');

% Filepaths to CSVs that were exported from ROS bags
base_path = 'C:\Users\riojaske\Documents\magsteer\PreopPlanning\data\phantom';
filepaths_p_manual.force = fullfile(base_path, 'manual\phantom_manual_ea_trial1_force.csv');
filepaths_p_manual.force2 = fullfile(base_path, 'manual\phantom_manual_ea_trial2_force.csv');
filepaths_p_manual.force3 = fullfile(base_path, 'manual\phantom_manual_ea_trial3_force.csv');
filepaths_p_manual.force4 = fullfile(base_path, 'manual\phantom_manual_ea2_trial4_force.csv');

filepaths_p_unguided.force = fullfile(base_path, 'UG-EA\phantom_ug_ea_trial1_force.csv');
filepaths_p_unguided.smaract = fullfile(base_path, 'UG-EA\phantom_ug_ea_trial1_smaract.csv');
filepaths_p_unguided2.force = fullfile(base_path, 'UG-EA\phantom_ug_ea_trial2_force.csv');
filepaths_p_unguided2.smaract = fullfile(base_path, 'UG-EA\phantom_ug_ea_trial2_smaract.csv');
% filepaths_mag.smaract   = fullfile(base_path, 'MagSteering_2019-07-22_trial4-mag_smaract.csv');
% filepaths_mag.force     = fullfile(base_path, 'MagSteering_2019-07-22_trial4-mag_force.csv');
% filepaths_nomag.smaract = fullfile(base_path, 'MagSteering_2019-07-22_trial4-nomag_smaract.csv');
% filepaths_nomag.force   = fullfile(base_path, 'MagSteering_2019-07-22_trial4-nomag_force.csv');

% Create data objects
data_manual = Nano17Data(filepaths_p_manual.force);
data_manual2 = Nano17Data(filepaths_p_manual.force2);
data_manual3 = Nano17Data(filepaths_p_manual.force3);
data_manual4 = Nano17Data(filepaths_p_manual.force4);
data_ug_ea = MagneticGuidanceData(filepaths_p_unguided);
data_ug2_ea = MagneticGuidanceData(filepaths_p_unguided2);

% set smoothing span (proportion of data points)
data_ug_ea.setSmoothSpan(0.1);
data_manual.setSmoothSpan(0.1);
%% plot
alpha = 1; % reduce transparency of unguided plot lines

% Figure: force components
h_forcecompare(1) = figure(1); clf;

subplot(4,2,1); grid on; hold on;
title('Manual');
ylim([-400,100]);
% manual, no smoothing
lh_manual(1) = plot(data_manual.time-data_manual.time(1),   data_manual.Fx, 'Color', [colors('jungle green'),  0.3*alpha]);
lh_manual(2) = plot(data_manual.time-data_manual.time(1),   data_manual.Fy, 'Color', [colors('red'), 0.3*alpha]);
lh_manual(3) = plot(data_manual.time-data_manual.time(1),   data_manual.Fz, 'Color', [colors('dark pastel purple'),   0.3*alpha]);

% manual, smoothed
lh_manual(4) = plot(data_manual.time-data_manual.time(1),   data_manual.Fx_smooth, 'Color', colors('jungle green'));
lh_manual(5) = plot(data_manual.time-data_manual.time(1),   data_manual.Fy_smooth, 'Color', colors('red'));
lh_manual(6) = plot(data_manual.time-data_manual.time(1),   data_manual.Fz_smooth, 'Color', colors('dark pastel purple'));

subplot(4,2,3); grid on; hold on;
title('Manual');
ylim([-400,100]);

% manual, no smoothing
lh_manual2(1) = plot(data_manual2.time-data_manual2.time(1),   data_manual2.Fx, 'Color', [colors('jungle green'),  0.3*alpha]);
lh_manual2(2) = plot(data_manual2.time-data_manual2.time(1),   data_manual2.Fy, 'Color', [colors('red'), 0.3*alpha]);
lh_manual2(3) = plot(data_manual2.time-data_manual2.time(1),   data_manual2.Fz, 'Color', [colors('dark pastel purple'),   0.3*alpha]);

% manual, smoothed
lh_manual2(4) = plot(data_manual2.time-data_manual2.time(1),   data_manual2.Fx_smooth, 'Color', colors('jungle green'));
lh_manual2(5) = plot(data_manual2.time-data_manual2.time(1),   data_manual2.Fy_smooth, 'Color', colors('red'));
lh_manual2(6) = plot(data_manual2.time-data_manual2.time(1),   data_manual2.Fz_smooth, 'Color', colors('dark pastel purple'));

subplot(4,2,5); grid on; hold on;
title('Manual');
ylim([-400,100]);

% manual, no smoothing
lh_manual3(1) = plot(data_manual3.time-data_manual3.time(1),   data_manual3.Fx, 'Color', [colors('jungle green'),  0.3*alpha]);
lh_manual3(2) = plot(data_manual3.time-data_manual3.time(1),   data_manual3.Fy, 'Color', [colors('red'), 0.3*alpha]);
lh_manual3(3) = plot(data_manual3.time-data_manual3.time(1),   data_manual3.Fz, 'Color', [colors('dark pastel purple'),   0.3*alpha]);

% manual, smoothed
lh_manual3(4) = plot(data_manual3.time-data_manual3.time(1),   data_manual3.Fx_smooth, 'Color', colors('jungle green'));
lh_manual3(5) = plot(data_manual3.time-data_manual3.time(1),   data_manual3.Fy_smooth, 'Color', colors('red'));
lh_manual3(6) = plot(data_manual3.time-data_manual3.time(1),   data_manual3.Fz_smooth, 'Color', colors('dark pastel purple'));

subplot(4,2,7); grid on; hold on;
title('Manual');
ylim([-400,100]);

% manual, no smoothing
lh_manual4(1) = plot(data_manual4.time-data_manual4.time(1),   data_manual4.Fx, 'Color', [colors('jungle green'),  0.3*alpha]);
lh_manual4(2) = plot(data_manual4.time-data_manual4.time(1),   data_manual4.Fy, 'Color', [colors('red'), 0.3*alpha]);
lh_manual4(3) = plot(data_manual4.time-data_manual4.time(1),   data_manual4.Fz, 'Color', [colors('dark pastel purple'),   0.3*alpha]);

% manual, smoothed
lh_manual4(4) = plot(data_manual4.time-data_manual4.time(1),   data_manual4.Fx_smooth, 'Color', colors('jungle green'));
lh_manual4(5) = plot(data_manual4.time-data_manual4.time(1),   data_manual4.Fy_smooth, 'Color', colors('red'));
lh_manual4(6) = plot(data_manual4.time-data_manual4.time(1),   data_manual4.Fz_smooth, 'Color', colors('dark pastel purple'));

subplot(4,2,2); grid on; hold on;

% unguided, no smoothing
lh_ug_ea(1) = plot(data_ug_ea.depth_insertion, data_ug_ea.Fx, 'Color', [colors('light green'),  0.3*alpha]);
lh_ug_ea(2) = plot(data_ug_ea.depth_insertion, data_ug_ea.Fy, 'Color', [colors('burnt orange'), 0.3*alpha]);
lh_ug_ea(3) = plot(data_ug_ea.depth_insertion, data_ug_ea.Fz, 'Color', [colors('light blue'),   0.3*alpha]);

% unguided, smoothed
lh_ug_ea(4) = plot(data_ug_ea.depth_insertion, data_ug_ea.Fx_smooth, 'Color', [colors('light green'), alpha]);
lh_ug_ea(5) = plot(data_ug_ea.depth_insertion, data_ug_ea.Fy_smooth, 'Color', [colors('burnt orange'), alpha]);
lh_ug_ea(6) = plot(data_ug_ea.depth_insertion, data_ug_ea.Fz_smooth, 'Color', [colors('light blue'), alpha]);

xlim([data_ug_ea.depth_insertion(1)-0.5, data_ug_ea.depth_insertion(end)+0.5]);
ylim([-400,100]);
title('UG-EA')

subplot(3,2,4); grid on; hold on;
% unguided, no smoothing
lh_ug2_ea(1) = plot(data_ug2_ea.depth_insertion, data_ug2_ea.Fx, 'Color', [colors('light green'),  0.3*alpha]);
lh_ug2_ea(2) = plot(data_ug2_ea.depth_insertion, data_ug2_ea.Fy, 'Color', [colors('burnt orange'), 0.3*alpha]);
lh_ug2_ea(3) = plot(data_ug2_ea.depth_insertion, data_ug2_ea.Fz, 'Color', [colors('light blue'),   0.3*alpha]);

% unguided, smoothed
lh_ug2_ea(4) = plot(data_ug2_ea.depth_insertion, data_ug2_ea.Fx_smooth, 'Color', [colors('light green'), alpha]);
lh_ug2_ea(5) = plot(data_ug2_ea.depth_insertion, data_ug2_ea.Fy_smooth, 'Color', [colors('burnt orange'), alpha]);
lh_ug2_ea(6) = plot(data_ug2_ea.depth_insertion, data_ug2_ea.Fz_smooth, 'Color', [colors('light blue'), alpha]);

set(lh_manual(1:3), 'LineWidth', 0.8)
set(lh_manual(4:6), 'LineWidth', 1.2)
set(lh_manual2(1:3), 'LineWidth', 0.8)
set(lh_manual2(4:6), 'LineWidth', 1.2)
set(lh_manual3(1:3), 'LineWidth', 0.8)
set(lh_manual3(4:6), 'LineWidth', 1.2)
set(lh_manual4(1:3), 'LineWidth', 0.8)
set(lh_manual4(4:6), 'LineWidth', 1.2)
set(lh_ug_ea(1:3), 'LineWidth', 0.8)
set(lh_ug_ea(4:6), 'LineWidth', 1)
set(lh_ug2_ea(1:3), 'LineWidth', 0.8)
set(lh_ug2_ea(4:6), 'LineWidth', 1)

% set(gca,'FontSize',15)
xlim([data_ug_ea.depth_insertion(1)-0.5, data_ug_ea.depth_insertion(end)+0.5]);
ylim([-400,100]);
title('UG-EA')
xlabel('Insertion Depth [mm]')
ylabel('Force [mN]')
% [~,icons,~,~] = legend([lh_ug_ea(1), lh_manual(1), lh_ug_ea(2), lh_manual(2), lh_ug_ea(3), lh_manual(3)], ...
%         'Fx','Fx_G','Fy','Fy_G','Fz','Fz_G', 'Location','southwest');
% set(icons(:),'LineWidth',2);

% Figure: Force Magnitude
% h_forcecompare(2) = figure(2); clf(2);
% grid on; hold on;
% lh2_ug_ea(1) = plot(data_ug_ea.depth_insertion, data_ug_ea.Fmag, 'Color', [colors('cocoa brown'),  0.3*alpha], 'LineWidth',1);
% lh2_ug_ea(2) = plot(data_ug_ea.depth_insertion, data_ug_ea.Fmag_smooth, 'Color', [colors('cocoa brown'),  alpha], 'LineWidth',2);
% lh2_manual(1) = plot(data_manual.time-data_manual.time(1), data_manual.Fmag, 'Color', [colors('dark lavender'), 0.3*alpha], 'LineWidth', 1);
% lh2_manual(2) = plot(data_manual.time-data_manual.time(1), data_manual.Fmag_smooth, 'Color', colors('dark lavender'), 'LineWidth', 2);
% % nomag_Fmag_interp = interp1(data_ug_ea.depth_insertion, data_ug_ea.Fmag_smooth, data_manual.depth_insertion);
% % lh2_diff = plot(data_manual.time-data_manual.time(1), (data_manual.Fmag_smooth - nomag_Fmag_interp), 'Color', colors('blue'), 'LineWidth', 2);
% set(gca,'FontSize',15)
% xlim([data_ug_ea.depth_insertion(1)-0.5, data_ug_ea.depth_insertion(end)+0.5])
% ylim([-inf,400]);
% title('Guided vs Unguided Insertion')
% xlabel('Insertion Depth [mm]')
% ylabel('||Force|| [mN]')

% [~,icons,~,~] = legend([lh2_ug_ea(2), lh2_manual(2), lh2_diff], 'Unguided', 'Guided', 'Difference', 'Location','northwest');
% [~,icons,~,~] = legend([lh2_ug_ea(2), lh2_manual(2)], 'Unguided', 'Guided', 'Location','northwest');set(icons(:),'LineWidth',2);
% savefig(h_forcecompare, 'CompareForces.fig')