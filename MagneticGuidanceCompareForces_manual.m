%% Data Analysis for Manual Insertion Trials
% Trevor Bruns and Katy Riojas
% Last Updated: 9/5/19

clear all; clc;

%% INPUTS
addpath('functions','classDef','data');
load('avg_cal_slope.mat');
smooth_span = 0.01;

% Specify filepaths to CSVs that were exported from ROS bags
base_path = 'data';

filepaths_p_manual.force = fullfile(base_path, 'phantom\manual\phantom_manual_ea_trial1_force.csv');
filepaths_p_manual.force2 = fullfile(base_path, 'phantom\manual\phantom_manual_ea_trial2_force.csv');
filepaths_p_manual.force3 = fullfile(base_path, 'phantom\manual\phantom_manual_ea_trial3_force.csv');
filepaths_p_manual.force4 = fullfile(base_path, 'phantom\manual\phantom_manual_ea2_trial4_force.csv');

filepaths_c_manual.force = fullfile(base_path, 'cadaver\manual\cadaver_manual_ea3_trial1_force.csv');
filepaths_c_manual.force2 = fullfile(base_path, 'cadaver\manual\cadaver_manual_ea3_trial2_force.csv');
filepaths_c_manual.force3 = fullfile(base_path, 'cadaver\manual\cadaver_manual_ea3_trial3_force.csv');

% Create data objects
data_pman1 = Nano17Data(filepaths_p_manual.force,avg_cal_slope);
data_pman2 = Nano17Data(filepaths_p_manual.force2,avg_cal_slope);
data_pman3 = Nano17Data(filepaths_p_manual.force3,avg_cal_slope);
data_pman4 = Nano17Data(filepaths_p_manual.force4,avg_cal_slope);

data_cman1 = Nano17Data(filepaths_c_manual.force,avg_cal_slope);
data_cman2 = Nano17Data(filepaths_c_manual.force2,avg_cal_slope);
data_cman3 = Nano17Data(filepaths_c_manual.force3,avg_cal_slope);

% Smooth Data -- note smoothing with time and force data now instead of
% depth and force data
data_pman1 = data_pman1.setSmoothSpan(smooth_span);
data_pman2 = data_pman2.setSmoothSpan(smooth_span);
data_pman3 = data_pman3.setSmoothSpan(smooth_span);
data_pman4 = data_pman4.setSmoothSpan(smooth_span);

data_cman1 = data_cman1.setSmoothSpan(smooth_span);
data_cman2 = data_cman2.setSmoothSpan(smooth_span);
data_cman3 = data_cman3.setSmoothSpan(smooth_span);
%% 
% Notes on Manual Insertion Trim Times From Cadaver Manual 
% Trial 1
tend_vid_peak = 69; % [s]
tend_F_peak = 81.8599; %[s]
t1_endinsertion_vid = 54; % [s]
diff_vid = tend_vid_peak - t1_endinsertion_vid; %[s]
t1start_F = 8.7; % [s]

% Trial 2
t2_vid_peak = 5.17; % first time force goes back to zero
t2_startinsertion_vid = 12.2; % first time when inserting
t2_endinsertion_vid = 83;
t2_start_F = 1;

% Trial 3
t3_peak_v = 3.68; % [s]
t3_startinsertion_v = 11.68; %[s]
t3start_F = 10.7;
tend_vid_end = 96.8; % [s]

%% Find trim points for time vectors
% visually determine stop and end points 
% Note - (labadie did large force spike at the beginning to indicate 
% about to begin start so we need to filter that out)
% Also trim end 
p1_start = find(data_pman1.time_vector>20,1);
p1_end= find(data_pman1.time_vector>90,1);
p1_vec = p1_start:p1_end;
p1_tvec = data_pman1.time_vector(p1_vec)-data_pman1.time_vector(p1_vec(1));

p2_start = find(data_pman2.time_vector>8.5,1);
p2_end= find(data_pman2.time_vector>146,1);
p2_vec = p2_start:p2_end;
p2_tvec = data_pman2.time_vector(p2_vec)-data_pman2.time_vector(p2_vec(1));

p3_start = find(data_pman3.time_vector>10,1);
p3_end= find(data_pman3.time_vector>88,1);
p3_vec = p3_start:p3_end;
p3_tvec = data_pman3.time_vector(p3_vec)-data_pman3.time_vector(p3_vec(1));

p4_start = find(data_pman4.time_vector>19,1);
p4_end= find(data_pman4.time_vector>163,1);
p4_vec = p4_start:p4_end;
p4_tvec = data_pman4.time_vector(p4_vec)-data_pman4.time_vector(p4_vec(1));

c1_start = find(data_cman1.time_vector>8,1);
c1_end= find(data_cman1.time_vector>78,1);
c1_vec = c1_start:c1_end;
c1_tvec = data_cman1.time_vector(c1_vec)-data_cman1.time_vector(c1_vec(1));

c2_start = find(data_cman2.time_vector>10,1);
c2_end= find(data_cman2.time_vector>90,1);
c2_vec = c2_start:c2_end;
c2_tvec = data_cman2.time_vector(c2_vec)-data_cman2.time_vector(c2_vec(1));

c3_start = find(data_cman3.time_vector>7,1);
c3_end= find(data_cman3.time_vector>103,1);
c3_vec = c3_start:c3_end;
c3_tvec = data_cman3.time_vector(c3_vec)-data_cman3.time_vector(c3_vec(1));

%% Plot 1: Raw XYZ Insertion Force Data vs. Time
alpha = 1; % line opacity
xyzColor = [1,0,0;0.1961,0.6275,0.2745;0,0,1];
cMat = distinguishable_colors(7);

figure(1); clf(1);
subplot(2,4,1); grid on; hold on; 
xlabel('Time (s)'); ylabel({'Phantom Trials';'Force [mN]'}); 
title('Trial 1'); 
xlim([0 70]); ylim([-250,60]);

plot(p1_tvec,data_pman1.Fx(p1_vec),'Color', xyzColor(1,:));
plot(p1_tvec,data_pman1.Fy(p1_vec),'Color', xyzColor(2,:));
plot(p1_tvec,data_pman1.Fz(p1_vec),'Color', xyzColor(3,:));

subplot(2,4,2); grid on; hold on;
title('Trial 2'); 
xlabel('Time (s)'); 
xlim([0 137]); ylim([-250,60]);
plot(p2_tvec,data_pman2.Fx(p2_vec),'Color', xyzColor(1,:));
plot(p2_tvec,data_pman2.Fy(p2_vec),'Color', xyzColor(2,:));
plot(p2_tvec,data_pman2.Fz(p2_vec),'Color', xyzColor(3,:));

subplot(2,4,3); grid on; hold on;
title('Trial 3');
xlabel('Time (s)'); 
xlim([0 77]); ylim([-250,60]);
plot(p3_tvec,data_pman3.Fx(p3_vec),'Color', xyzColor(1,:));
plot(p3_tvec,data_pman3.Fy(p3_vec),'Color', xyzColor(2,:));
plot(p3_tvec,data_pman3.Fz(p3_vec),'Color', xyzColor(3,:));

subplot(2,4,4); grid on; hold on;
title('Trial 4'); 
xlabel('Time (s)');
xlim([0 144]); ylim([-250,60]);
plot(p4_tvec,data_pman4.Fx(p4_vec),'Color', xyzColor(1,:));
plot(p4_tvec,data_pman4.Fy(p4_vec),'Color', xyzColor(2,:));
plot(p4_tvec,data_pman4.Fz(p4_vec),'Color', xyzColor(3,:));

subplot(2,3,4); grid on; hold on;
title('Trial 1');
xlim([0 67]); ylim([-250,45]);
xlabel('Time (s)'); ylabel({'Cadaver Trials';'Force [mN]'}); 
plot(c1_tvec,data_cman1.Fx(c1_vec),'Color', xyzColor(1,:));
plot(c1_tvec,data_cman1.Fy(c1_vec),'Color', xyzColor(2,:));
plot(c1_tvec,data_cman1.Fz(c1_vec),'Color', xyzColor(3,:));

subplot(2,3,5); grid on; hold on;
xlabel('Time (s)');
xlim([0 80]); ylim([-250,45]);
title('Trial 2');
plot(c2_tvec,data_cman2.Fx(c2_vec),'Color', xyzColor(1,:));
plot(c2_tvec,data_cman2.Fy(c2_vec),'Color', xyzColor(2,:));
plot(c2_tvec,data_cman2.Fz(c2_vec),'Color', xyzColor(3,:));

subplot(2,3,6); grid on; hold on;
xlabel('Time (s)');
xlim([0 95]); ylim([-250,45]);
title('Trial 3');
plot(c3_tvec,data_cman3.Fx(c3_vec),'Color', xyzColor(1,:));
plot(c3_tvec,data_cman3.Fy(c3_vec),'Color', xyzColor(2,:));
plot(c3_tvec,data_cman3.Fz(c3_vec),'Color', xyzColor(3,:));
legend('Fx','Fy','Fz');

%% Plot 2: Compare Magnitude across trials
figure(2); clf(2);
subplot(1,2,1); grid on; hold on;
xlabel('Time (s)'); ylabel('Force [mN]');
title('Phantom Trials: Fmag vs. Time');

plot(data_pman1.time_vector, data_pman1.Fmag, 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h1 = plot(data_pman1.time_vector, data_pman1.Fmag_smooth, 'Color', cMat(1,:), 'LineWidth',1);

plot(data_pman2.time_vector, data_pman2.Fmag, 'Color', [cMat(2,:),  0.3*alpha], 'LineWidth',0.1);
h2 = plot(data_pman2.time_vector, data_pman2.Fmag_smooth, 'Color', cMat(2,:), 'LineWidth',1);

plot(data_pman3.time_vector, data_pman3.Fmag, 'Color', [cMat(2,:),  0.3*alpha], 'LineWidth',0.1);
h3 = plot(data_pman3.time_vector, data_pman3.Fmag_smooth, 'Color', cMat(3,:), 'LineWidth',1);

plot(data_pman4.time_vector, data_pman4.Fmag, 'Color', [cMat(4,:),  0.3*alpha], 'LineWidth',0.1);
h4 = plot(data_pman4.time_vector, data_pman4.Fmag_smooth, 'Color', cMat(4,:), 'LineWidth',1);

legend([h1,h2,h3,h4],{'Trial 1','Trial 2','Trial 3','Trial 4'});

subplot(1,2,2); grid on; hold on;
xlabel('Time (s)'); ylabel('Force [mN]');
title('Cadaver Trials: Fmag vs. Time');

plot(data_cman1.time_vector, data_cman1.Fmag, 'Color', [cMat(5,:),  0.3*alpha], 'LineWidth',0.1);
h5 = plot(data_cman1.time_vector, data_cman1.Fmag_smooth, 'Color', cMat(5,:), 'LineWidth',1);

plot(data_cman2.time_vector, data_cman2.Fmag, 'Color', [cMat(6,:),  0.3*alpha], 'LineWidth',0.1);
h6 = plot(data_cman2.time_vector, data_cman2.Fmag_smooth, 'Color', cMat(6,:), 'LineWidth',1);

plot(data_cman3.time_vector, data_cman3.Fmag, 'Color', [cMat(7,:),  0.3*alpha], 'LineWidth',0.1);
h7= plot(data_cman3.time_vector, data_cman3.Fmag_smooth, 'Color', cMat(7,:), 'LineWidth',1);

legend([h5,h6,h7],{'Trial 1','Trial 2','Trial 3'});

%% Plot 3: Compare Magnitude of Trial 3 Phantom and Trial 2 Cadaver
figure(3); clf(3); 
subplot(1,2,1); grid on; hold on;
xlabel('Time (s)'); ylabel('Force Magnitude (mN)');
ylim([0 250]);
title('(a)');
praw = plot(p2_tvec,data_pman2.Fmag(p2_vec), 'Color', ['r',  0.3*alpha], 'LineWidth',0.1);
psmooth = plot(p2_tvec,data_pman2.Fmag_smooth(p2_vec), 'Color', 'r', 'LineWidth',1);
legend([praw,psmooth],'Raw','Smoothed');
    
subplot(1,2,2); grid on; hold on;
xlabel('Time (s)'); ylabel('Force Magnitude (mN)');
ylim([0 250]);
title('(b)');
craw = plot(c2_tvec,data_cman2.Fmag(c2_vec), 'Color', ['b',  0.3*alpha], 'LineWidth',0.1);
csmooth = plot(c2_tvec,data_cman2.Fmag_smooth(c2_vec), 'Color', 'b', 'LineWidth',1);
legend([craw,csmooth],'Raw','Smoothed');

%% Save Raw Force Magnitudes
pmanFmag = [data_pman1.Fmag;data_pman2.Fmag;data_pman3.Fmag;data_pman4.Fmag];
cmanFmag = [data_cman1.Fmag;data_cman2.Fmag;data_cman3.Fmag];

save('data\phantom\pmanFmag.mat','pmanFmag');
save('data\cadaver\cmanFmag.mat','cmanFmag');

%% Plot 4: Force magnitude vs. Time for cropping
LL = 5;
pman1_peak_idx = find(data_pman1.Fmag>50,1);
tPP_pman1_idx = find(data_pman1.Fmag(pman1_peak_idx:end)<LL,1);
tPP_pman1 = data_pman1.time_vector(tPP_pman1_idx+pman1_peak_idx);

pman2_peak_idx = find(data_pman2.Fmag>135.3,1);
tPP_pman2_idx = find(data_pman2.Fmag(pman2_peak_idx:end)<LL,1);
tPP_pman2 = data_pman2.time_vector(tPP_pman2_idx+pman2_peak_idx);

pman3_peak_idx = find(data_pman3.Fmag>50,1);
tPP_pman3_idx = find(data_pman3.Fmag(pman3_peak_idx:end)<LL,1);
tPP_pman3 = data_pman3.time_vector(tPP_pman3_idx+pman3_peak_idx);

pman4_peak_idx = find(data_pman4.Fmag>50,1);
tPP_pman4_idx = find(data_pman4.Fmag(pman4_peak_idx:end)<LL,1);
tPP_pman4 = data_pman4.time_vector(tPP_pman4_idx+pman4_peak_idx);

cman1_peak_idx = find(data_cman1.Fmag>50,1);
tPP_cman1_idx = find(data_cman1.Fmag(cman1_peak_idx:end)<LL,1);
tPP_cman1 = data_cman1.time_vector(tPP_cman1_idx+cman1_peak_idx);

cman2_peak_idx = find(data_cman2.Fmag>50,1);
tPP_cman2_idx = find(data_cman2.Fmag(cman2_peak_idx:end)<LL,1);
tPP_cman2 = data_cman2.time_vector(tPP_cman2_idx+cman2_peak_idx);

cman3_peak_idx = find(data_cman3.Fmag>50,1);
tPP_cman3_idx = find(data_cman3.Fmag(cman3_peak_idx:end)<LL,1);
tPP_cman3 = data_cman3.time_vector(tPP_cman3_idx+cman3_peak_idx);

figure(4); clf(4); 

subplot(2,4,1); grid on; hold on; xlabel('Time (s)'); ylabel ('||F|| (mN)');
plot(data_pman1.time_vector,data_pman1.Fmag,'Color','b');
subplot(2,4,2); grid on; hold on; 
plot(data_pman2.time_vector,data_pman2.Fmag,'Color','b');
subplot(2,4,3); grid on; hold on; 
plot(data_pman3.time_vector,data_pman3.Fmag,'Color','b');
subplot(2,4,4); grid on; hold on; 
plot(data_pman4.time_vector,data_pman4.Fmag,'Color','b');

subplot(2,3,4); grid on; hold on; 
plot(data_cman1.time_vector,data_cman1.Fmag,'Color','b');
subplot(2,3,5); grid on; hold on; 
plot(data_cman2.time_vector,data_cman2.Fmag,'Color','b');
subplot(2,3,6); grid on; hold on; 
plot(data_cman3.time_vector,data_cman3.Fmag,'Color','b');


