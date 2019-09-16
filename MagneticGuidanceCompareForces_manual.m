%% Data Analysis for Manual Insertion Trials
% Trevor Bruns and Katy Riojas
% Last Updated: 9/8/19

% Function "trimManualTrials" trims according to video footage

%% Import Force and Smaract data
clear all; clc;

addpath('functions','classDef','data');
% load force sensor calibration slopes
load('avg_cal_slopes.mat');

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

% save('ManualPData.mat','data_pman1','data_pman2','data_pman3','data_pman4');
%% Trim Trials
[p1_tvec,p1_vec,p2_tvec,p2_vec,p3_tvec,p3_vec,...
    p4_tvec,p4_vec,c1_tvec,c1_vec,c2_tvec,c2_vec,c3_tvec,c3_vec,...
    releaseTimes] = trimManualTrials(data_pman1,data_pman2,data_pman3,...
    data_pman4,data_cman1,data_cman2,data_cman3);

%% Plotting Variables
xyzColor = [1,0,0;0.1961,0.6275,0.2745;0,0,1];
alpha = 1; % reduce transparency of unguided plot lines
cMat = distinguishable_colors(4);
LW1 = 2;
% Calc max y phantom
maxpY = max([data_pman1.Fmag;data_pman2.Fmag;...
    data_pman3.Fmag;data_pman4.Fmag]);
% Calc max y cadaver
maxcY = max([data_cman1.Fmag;data_cman2.Fmag;data_cman3.Fmag]);
maxcY_trimmed = max([data_cman1.Fmag(c1_vec);data_cman2.Fmag(c2_vec);data_cman3.Fmag(c3_vec)]);
plotXYZData = 0;

%% Plot 1: Force magnitude vs. Time for cropping
figure(1); clf(1); 
sgtitle('Raw Fmag Data vs. Time all Trials');
subplot(2,4,1); grid on; hold on; xlabel('Time (s)'); 
ylabel ('Phantom ||F|| (mN)');
ylim([0,maxpY]);
plot(data_pman1.time,data_pman1.Fmag,'Color','b');

line([p1_tvec(1),p1_tvec(1)],[0,max(data_pman1.Fmag)],'Color','k','LineStyle','--');
line([p1_tvec(end),p1_tvec(end)],[0,max(data_pman1.Fmag)],'Color','k','LineStyle','--');
line([releaseTimes(1),releaseTimes(1)],[0,maxpY],'Color','g','LineStyle','--');
scatter(p1_tvec(1),data_pman1.Fmag(p1_vec(1)),10,'filled','r');
scatter(p1_tvec(end),data_pman2.Fmag(p1_vec(end)),10,'filled','r');

subplot(2,4,2); grid on; hold on; xlabel('Time (s)'); 
ylim([0,maxpY]);
plot(data_pman2.time,data_pman2.Fmag,'Color','b');
line([p2_tvec(1),p2_tvec(1)],[0,max(data_pman2.Fmag)],'Color','k','LineStyle','--');
line([p2_tvec(end),p2_tvec(end)],[0,max(data_pman2.Fmag)],'Color','k','LineStyle','--');
line([releaseTimes(2),releaseTimes(2)],[0,maxpY],'Color','g','LineStyle','--');
scatter(p2_tvec(1),data_pman2.Fmag(p2_vec(1)),10,'filled','r');
scatter(p2_tvec(end),data_pman2.Fmag(p2_vec(end)),10,'filled','r');

subplot(2,4,3); grid on; hold on; xlabel('Time (s)'); 
ylim([0,maxpY]);
plot(data_pman3.time,data_pman3.Fmag,'Color','b');
line([p3_tvec(1),p3_tvec(1)],[0,max(data_pman3.Fmag)],'Color','k','LineStyle','--');
line([p3_tvec(end),p3_tvec(end)],[0,max(data_pman3.Fmag)],'Color','k','LineStyle','--');
line([releaseTimes(3),releaseTimes(3)],[0,maxpY],'Color','g','LineStyle','--');
scatter(p3_tvec(1),data_pman3.Fmag(p3_vec(1)),10,'filled','r');
scatter(p3_tvec(end),data_pman3.Fmag(p3_vec(end)),10,'filled','r');

subplot(2,4,4); grid on; hold on; xlabel('Time (s)'); 
ylim([0,maxpY]);
plot(data_pman4.time,data_pman4.Fmag,'Color','b');
line([p4_tvec(1),p4_tvec(1)],[0,max(data_pman4.Fmag)],'Color','k','LineStyle','--');
line([p4_tvec(end),p4_tvec(end)],[0,max(data_pman4.Fmag)],'Color','k','LineStyle','--');
line([releaseTimes(4),releaseTimes(4)],[0,maxpY],'Color','g','LineStyle','--');
scatter(p4_tvec(1),data_pman4.Fmag(p4_vec(1)),10,'filled','r');
scatter(p4_tvec(end),data_pman4.Fmag(p4_vec(end)),10,'filled','r');

subplot(2,3,4); grid on; hold on; ylabel('Cadaver ||Force|| [mN]');
ylim([0,maxcY]);
plot(data_cman1.time,data_cman1.Fmag,'Color','b'); xlabel('Time (s)'); 
line([c1_tvec(1),c1_tvec(1)],[0,max(data_cman1.Fmag)],'Color','k','LineStyle','--');
line([c1_tvec(end),c1_tvec(end)],[0,max(data_cman1.Fmag)],'Color','k','LineStyle','--');
% line([releaseTimes(5),releaseTimes(5)],[0,maxpY],'Color','g','LineStyle','--');
scatter(c1_tvec(1),data_cman1.Fmag(c1_vec(1)),10,'filled','r');
scatter(c1_tvec(end),data_cman1.Fmag(c1_vec(end)),10,'filled','r');

subplot(2,3,5); grid on; hold on;
ylim([0,maxcY]);
plot(data_cman2.time,data_cman2.Fmag,'Color','b');
line([c2_tvec(1),c2_tvec(1)],[0,max(data_cman2.Fmag)],'Color','k','LineStyle','--');
line([c2_tvec(end),c2_tvec(end)],[0,max(data_cman2.Fmag)],'Color','k','LineStyle','--');
scatter(c2_tvec(1),data_cman2.Fmag(c2_vec(1)),10,'filled','r');
scatter(c2_tvec(end),data_cman2.Fmag(c2_vec(end)),10,'filled','r');

subplot(2,3,6); grid on; hold on; xlabel('Time (s)'); 
ylim([0,maxcY]);
plot(data_cman3.time,data_cman3.Fmag,'Color','b');
line([c3_tvec(1),c3_tvec(1)],[0,max(data_cman3.Fmag)],'Color','k','LineStyle','--');
line([c3_tvec(end),c3_tvec(end)],[0,max(data_cman3.Fmag)],'Color','k','LineStyle','--');
scatter(c3_tvec(1),data_cman3.Fmag(c3_vec(1)),10,'filled','r');
scatter(c3_tvec(end),data_cman3.Fmag(c3_vec(end)),10,'filled','r');

%% Plot 2: Compare Magnitude across trials with smoothing and trimming
figure(2); clf(2); 
sgtitle('Cropped Fmag Data (Smooth and Raw) vs. Time all Trials');
subplot(2,4,1); grid on; hold on; 
xlabel('Time (s)'); ylabel ('Phantom ||F|| (mN)');
title('Trial 1');
plot(p1_tvec, data_pman1.Fmag(p1_vec), 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h1 = plot(p1_tvec, data_pman1.Fmag_smooth(p1_vec), 'Color', [cMat(1,:),  alpha], 'LineWidth',1);

subplot(2,4,2); grid on; hold on; 
xlabel('Time (s)'); title('Trial 2');
plot(p2_tvec, data_pman2.Fmag(p2_vec), 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h2 = plot(p2_tvec, data_pman2.Fmag_smooth(p2_vec), 'Color', [cMat(1,:),  alpha], 'LineWidth',1);

subplot(2,4,3); grid on; hold on; 
xlabel('Time (s)'); title('Trial 3');
plot(p3_tvec, data_pman3.Fmag(p3_vec), 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h3 = plot(p3_tvec, data_pman3.Fmag_smooth(p3_vec), 'Color', [cMat(1,:),  alpha], 'LineWidth',1);

subplot(2,4,4); grid on; hold on; 
xlabel('Time (s)'); title('Trial 4');
plot(p4_tvec, data_pman4.Fmag(p4_vec), 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h4 = plot(p4_tvec, data_pman4.Fmag_smooth(p4_vec), 'Color', [cMat(1,:),  alpha], 'LineWidth',1);

subplot(2,3,4); grid on; hold on;
xlabel('Time (s)'); ylabel('Cadaver ||Force|| [mN]');
title('Trial 1');
plot(c1_tvec, data_cman1.Fmag(c1_vec), 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h5 = plot(c1_tvec, data_cman1.Fmag_smooth(c1_vec), 'Color', [cMat(1,:),  alpha], 'LineWidth',1);

subplot(2,3,5); grid on; hold on;
xlabel('Time (s)'); title('Trial 2');
plot(c2_tvec, data_cman2.Fmag(c2_vec), 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h6 = plot(c2_tvec, data_cman2.Fmag_smooth(c2_vec), 'Color', [cMat(1,:),  alpha], 'LineWidth',1);

subplot(2,3,6); grid on; hold on;
xlabel('Time (s)'); title('Trial 3');
plot(c3_tvec, data_cman3.Fmag(c3_vec), 'Color', [cMat(1,:),  0.3*alpha], 'LineWidth',0.1);
h7= plot(c3_tvec, data_cman3.Fmag_smooth(c3_vec), 'Color', [cMat(1,:),  alpha], 'LineWidth',1);

%% Plot the Force Magnitude vs. Time for just the Cadaver
figure(3); clf(3);
grid on; hold on; 
xlabel('Insertion Time (s)'); ylabel('Force (mN)');
ylim([0,120]);
c1Tvec = c1_tvec - c1_tvec(1);
c2Tvec = c2_tvec - c2_tvec(1);
c3Tvec = c3_tvec - c3_tvec(1);

xlim([0,max([c1Tvec;c2Tvec;c3Tvec])]);
plot(c1Tvec,data_cman1.Fmag_smooth(c1_vec),'Color','b','LineWidth',1);
plot(c2Tvec,data_cman2.Fmag_smooth(c2_vec),'Color','r','LineWidth',1);
plot(c3Tvec,data_cman3.Fmag_smooth(c3_vec),'Color','g','LineWidth',1);
legend('Manual Cadaver Trial 1','Manual Cadaver Trial 2','Manual Cadaver Trial 3');

%% Save Raw Force Magnitudes
pmanFmag = [data_pman1.Fmag(p1_vec);...
            data_pman2.Fmag(p2_vec);...
            data_pman3.Fmag(p3_vec);...
            data_pman4.Fmag(p4_vec)];
cmanFmag = [data_cman1.Fmag(c1_vec);...
            data_cman2.Fmag(c2_vec);...
            data_cman3.Fmag(c3_vec)];

save('data\phantom\pmanFmag.mat','pmanFmag');
save('data\cadaver\cmanFmag.mat','cmanFmag');

%% Plot 3: XYZ Plots using trimmed data
if plotXYZData
    figure(4); clf(4);
    subplot(2,4,1); grid on; hold on; 
    xlabel('Time (s)'); ylabel({'Phantom Trials';'Force [mN]'}); 
    sgtitle('XYZ Data from All Manual Trials');
    title('Trial 1'); 
    plot(p1_tvec,data_pman1.Fx(p1_vec),'Color', xyzColor(1,:));
    plot(p1_tvec,data_pman1.Fy(p1_vec),'Color', xyzColor(2,:));
    plot(p1_tvec,data_pman1.Fz(p1_vec),'Color', xyzColor(3,:));

    subplot(2,4,2); grid on; hold on;
    title('Trial 2'); xlabel('Time (s)'); 
    plot(p2_tvec,data_pman2.Fx(p2_vec),'Color', xyzColor(1,:));
    plot(p2_tvec,data_pman2.Fy(p2_vec),'Color', xyzColor(2,:));
    plot(p2_tvec,data_pman2.Fz(p2_vec),'Color', xyzColor(3,:));

    subplot(2,4,3); grid on; hold on;
    title('Trial 3'); xlabel('Time (s)'); 
    plot(p3_tvec,data_pman3.Fx(p3_vec),'Color', xyzColor(1,:));
    plot(p3_tvec,data_pman3.Fy(p3_vec),'Color', xyzColor(2,:));
    plot(p3_tvec,data_pman3.Fz(p3_vec),'Color', xyzColor(3,:));

    subplot(2,4,4); grid on; hold on;
    title('Trial 4'); xlabel('Time (s)');
    plot(p4_tvec,data_pman4.Fx(p4_vec),'Color', xyzColor(1,:));
    plot(p4_tvec,data_pman4.Fy(p4_vec),'Color', xyzColor(2,:));
    plot(p4_tvec,data_pman4.Fz(p4_vec),'Color', xyzColor(3,:));

    subplot(2,3,4); grid on; hold on;
    title('Trial 1');
    xlabel('Time (s)'); ylabel({'Cadaver Trials';'Force [mN]'}); 
    plot(c1_tvec,data_cman1.Fx(c1_vec),'Color', xyzColor(1,:));
    plot(c1_tvec,data_cman1.Fy(c1_vec),'Color', xyzColor(2,:));
    plot(c1_tvec,data_cman1.Fz(c1_vec),'Color', xyzColor(3,:));

    subplot(2,3,5); grid on; hold on;
    xlabel('Time (s)');
    title('Trial 2');
    plot(c2_tvec,data_cman2.Fx(c2_vec),'Color', xyzColor(1,:));
    plot(c2_tvec,data_cman2.Fy(c2_vec),'Color', xyzColor(2,:));
    plot(c2_tvec,data_cman2.Fz(c2_vec),'Color', xyzColor(3,:));

    subplot(2,3,6); grid on; hold on;
    xlabel('Time (s)');
    title('Trial 3');
    plot(c3_tvec,data_cman3.Fx(c3_vec),'Color', xyzColor(1,:));
    plot(c3_tvec,data_cman3.Fy(c3_vec),'Color', xyzColor(2,:));
    plot(c3_tvec,data_cman3.Fz(c3_vec),'Color', xyzColor(3,:));
    legend('Fx','Fy','Fz');
end