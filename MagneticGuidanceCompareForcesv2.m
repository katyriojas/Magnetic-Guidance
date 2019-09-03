%% Import Force and Smaract data
clear all; close all; clc;

addpath('C:\Users\riojaske\Documents\magsteer\Magnetic-Guidance\functions');
addpath('C:\Users\riojaske\Documents\magsteer\Magnetic-Guidance\classDef');

% Filepaths to CSVs that were exported from ROS bags
base_path = 'C:\Users\riojaske\Documents\magsteer\Magnetic-Guidance\data';

filepaths_p_manual.force = fullfile(base_path, 'phantom\manual\phantom_manual_ea_trial1_force.csv');
filepaths_p_manual.force2 = fullfile(base_path, 'phantom\manual\phantom_manual_ea_trial2_force.csv');
filepaths_p_manual.force3 = fullfile(base_path, 'phantom\manual\phantom_manual_ea_trial3_force.csv');
filepaths_p_manual.force4 = fullfile(base_path, 'phantom\manual\phantom_manual_ea2_trial4_force.csv');

filepaths_c_manual.force = fullfile(base_path, 'cadaver\manual\cadaver_manual_ea3_trial1_force.csv');
filepaths_c_manual.force2 = fullfile(base_path, 'cadaver\manual\cadaver_manual_ea3_trial2_force.csv');
filepaths_c_manual.force3 = fullfile(base_path, 'cadaver\manual\cadaver_manual_ea3_trial3_force.csv');

% Create data objects
data_pman1 = Nano17Data(filepaths_p_manual.force);
data_pman2 = Nano17Data(filepaths_p_manual.force2);
data_pman3 = Nano17Data(filepaths_p_manual.force3);
data_pman4 = Nano17Data(filepaths_p_manual.force4);

data_cman1 = Nano17Data(filepaths_c_manual.force);
data_cman2 = Nano17Data(filepaths_c_manual.force2);
data_cman3 = Nano17Data(filepaths_c_manual.force3);

smooth_span = 0.01;

pman1_Fmag_smooth = smoothForces(data_pman1,smooth_span);
pman2_Fmag_smooth = smoothForces(data_pman2,smooth_span);
pman3_Fmag_smooth = smoothForces(data_pman3,smooth_span);
pman4_Fmag_smooth = smoothForces(data_pman4,smooth_span);

cman1_Fmag_smooth = smoothForces(data_cman1,smooth_span);
cman2_Fmag_smooth = smoothForces(data_cman2,smooth_span);
cman3_Fmag_smooth = smoothForces(data_cman3,smooth_span);

%% Plot the manual data
alpha = 1; % reduce transparency of unguided plot lines
cMat = distinguishable_colors(7);

figure(1); 
subplot(1,2,1); grid on; hold on;
xlabel('Time (s)'); ylabel('Force [mN]');
title('Manual Insertions in Phantoms');

% manual, no smoothing
plot(data_pman1.time-data_pman1.time(1),data_pman1.Fx, 'Color', [cMat(1,:), 0.3*alpha]);
plot(data_pman1.time-data_pman1.time(1),data_pman1.Fy, 'Color', [cMat(1,:), 0.3*alpha]);
plot(data_pman1.time-data_pman1.time(1),data_pman1.Fz, 'Color', [cMat(1,:), 0.3*alpha]);

plot(data_pman2.time-data_pman2.time(1),data_pman2.Fx, 'Color', [cMat(2,:), 0.3*alpha]);
plot(data_pman2.time-data_pman2.time(1),data_pman2.Fy, 'Color', [cMat(2,:), 0.3*alpha]);
plot(data_pman2.time-data_pman2.time(1),data_pman2.Fz, 'Color', [cMat(2,:), 0.3*alpha]);

plot(data_pman3.time-data_pman3.time(1),data_pman3.Fx, 'Color', [cMat(3,:), 0.3*alpha]);
plot(data_pman3.time-data_pman3.time(1),data_pman3.Fy, 'Color', [cMat(3,:), 0.3*alpha]);
plot(data_pman3.time-data_pman3.time(1),data_pman3.Fz, 'Color', [cMat(3,:), 0.3*alpha]);

plot(data_pman4.time-data_pman4.time(1),data_pman4.Fx, 'Color', [cMat(4,:), 0.3*alpha]);
plot(data_pman4.time-data_pman4.time(1),data_pman4.Fy, 'Color', [cMat(4,:), 0.3*alpha]);
plot(data_pman4.time-data_pman4.time(1),data_pman4.Fz, 'Color', [cMat(4,:), 0.3*alpha]);

subplot(1,2,2); grid on; hold on;
xlabel('Time (s)'); ylabel('Force [mN]');
title('Manual Insertions in a Cadaver');

plot(data_cman1.time-data_cman1.time(1),data_cman1.Fx, 'Color', [cMat(5,:), 0.3*alpha]);
plot(data_cman1.time-data_cman1.time(1),data_cman1.Fy, 'Color', [cMat(5,:), 0.3*alpha]);
plot(data_cman1.time-data_cman1.time(1),data_cman1.Fz, 'Color', [cMat(5,:), 0.3*alpha]);

plot(data_cman2.time-data_cman2.time(1),data_cman2.Fx, 'Color', [cMat(6,:), 0.3*alpha]);
plot(data_cman2.time-data_cman2.time(1),data_cman2.Fy, 'Color', [cMat(6,:), 0.3*alpha]);
plot(data_cman2.time-data_cman2.time(1),data_cman2.Fz, 'Color', [cMat(6,:), 0.3*alpha]);

plot(data_cman3.time-data_cman3.time(1),data_cman3.Fx, 'Color', [cMat(7,:), 0.3*alpha]);
plot(data_cman3.time-data_cman3.time(1),data_cman3.Fy, 'Color', [cMat(7,:), 0.3*alpha]);
plot(data_cman3.time-data_cman3.time(1),data_cman3.Fz, 'Color', [cMat(7,:), 0.3*alpha]);

figure(2); 
subplot(1,2,1); grid on; hold on;
xlabel('Time (s)'); ylabel('Force [mN]');
title('Magnitude of Manual Insertions in Phantoms');
% plot(10^-9*(data_pman1.time-data_pman1.time(1)), data_pman1.Fmag, 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h1 = plot(10^-9*(data_pman1.time-data_pman1.time(1)), pman1_Fmag_smooth, 'Color', [cMat(1,:),  alpha], 'LineWidth',1);

% plot(10^-9*(data_pman2.time-data_pman2.time(1)), data_pman2.Fmag, 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h2 = plot(10^-9*(data_pman2.time-data_pman2.time(1)), pman2_Fmag_smooth, 'Color', [cMat(2,:),  0.3*alpha], 'LineWidth',1);

% plot(10^-9*(data_pman3.time-data_pman3.time(1)), data_pman3.Fmag, 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h3 = plot(10^-9*(data_pman3.time-data_pman3.time(1)), pman3_Fmag_smooth, 'Color', [cMat(3,:),  0.3*alpha], 'LineWidth',1);

% plot(10^-9*(data_pman4.time-data_pman4.time(1)), data_pman4.Fmag, 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h4 = plot(10^-9*(data_pman4.time-data_pman4.time(1)), pman4_Fmag_smooth, 'Color', [cMat(4,:),  0.3*alpha], 'LineWidth',1);

legend([h1,h2,h3,h4],{'Trial 1','Trial 2','Trial 3','Trial 4'});

subplot(1,2,2); grid on; hold on;
xlabel('Time (s)'); ylabel('Force [mN]');
title('Magnitude of Manual Insertions in Cadavers');
% plot(10^-9*(data_cman1.time-data_cman1.time(1)), data_cman1.Fmag, 'Color', [cMat(5,:),  0.3*alpha], 'LineWidth',0.1);
h5 = plot(10^-9*(data_cman1.time-data_cman1.time(1)), cman1_Fmag_smooth, 'Color', [cMat(5,:),  0.3*alpha], 'LineWidth',1);

% plot(10^-9*(data_cman2.time-data_cman2.time(1)), data_cman2.Fmag, 'Color', [cMat(6,:),  0.3*alpha], 'LineWidth',0.1);
h6 = plot(10^-9*(data_cman2.time-data_cman2.time(1)), cman2_Fmag_smooth, 'Color', [cMat(6,:),  0.3*alpha], 'LineWidth',1);

% plot(10^-9*(data_cman3.time-data_cman3.time(1)), data_cman3.Fmag, 'Color', [cMat(7,:),  0.3*alpha], 'LineWidth',0.1);
h7= plot(10^-9*(data_cman3.time-data_cman3.time(1)), cman3_Fmag_smooth, 'Color', [cMat(7,:),  0.3*alpha], 'LineWidth',1);

legend([h5,h6,h7],{'Trial 1','Trial 2','Trial 3'});

xvec = linspace(0,100,1000);

Fmag_smooth_c = [interp1(10^-9*(data_cman1.time-data_cman1.time(1)),cman1_Fmag_smooth,xvec);...
                 interp1(10^-9*(data_cman2.time-data_cman2.time(1)),cman2_Fmag_smooth,xvec);...
                 interp1(10^-9*(data_cman3.time-data_cman3.time(1)),cman3_Fmag_smooth,xvec)];

Favg_manual_c = nanmean(Fmag_smooth_c,1);
std_manual_c = nanstd(Fmag_smooth_c);

figure(4); 
grid on; hold on; xlabel('Time (s)'); ylabel('Force (mN)');
title('Average Insertion Force in Cadaver vs. Time');
plot(xvec,Favg_manual_c,'Color',cMat(1,:));
fill([xvec fliplr(xvec)],[Favg_manual_c + std_manual_c, fliplr(Favg_manual_c-std_manual_c)],...
    'b','FaceAlpha',0.2,'LineStyle','none');
