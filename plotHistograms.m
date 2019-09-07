% Katy Riojas
% 9/6/19
clear all; clc;

% Plot a histogram of Force Percentages during insertions
% Pull in manual data
load('data\phantom\pmanFmag.mat'); % load pmanFmag
load('data\cadaver\cmanFmag.mat'); % load cmanFmag

load('data\phantom\p_ug_Fmag.mat','p_ug_Fmag'); % load p_ug_Fmag
load('data\phantom\p_g_Fmag.mat','p_g_Fmag'); % load p_g_Fmag

load('data\cadaver\c_ug_Fmag.mat','c_ug_Fmag'); % load c_ug_Fmag
load('data\cadaver\c_g_Fmag.mat','c_g_Fmag'); % load c_g_Fmag

% numbins = 50;
binWdth = 1;
xlimitMaxP = max([pmanFmag;p_ug_Fmag;p_g_Fmag]);
xFitP = linspace(0,xlimitMaxP,500);
xlimitMaxC = max([cmanFmag;c_ug_Fmag;c_g_Fmag]);
xFitC = linspace(0,xlimitMaxC,500);

% Fit gamma
pmanFmag_dist = fitdist(sort(pmanFmag),'gamma');
pmanFmag_dist.ParameterValues
pmanFmag_Y = gampdf(xFitP,pmanFmag_dist.ParameterValues(1),pmanFmag_dist.ParameterValues(2));

p_ug_Fmag_dist = fitdist(sort(p_ug_Fmag),'gamma');
p_ug_Fmag_dist.ParameterValues
p_ug_Fmag_Y = gampdf(xFitP,p_ug_Fmag_dist.ParameterValues(1),p_ug_Fmag_dist.ParameterValues(2));

p_g_Fmag_dist = fitdist(sort(p_g_Fmag),'gamma');
p_g_Fmag_dist.ParameterValues
p_g_Fmag_Y = gampdf(xFitP,p_g_Fmag_dist.ParameterValues(1),p_g_Fmag_dist.ParameterValues(2));

figure(4); clf(4);
subplot(1,2,1); grid on; hold on;
title('Phantom Insertions');
xlabel('||Force|| (mN)'); ylabel('Percentage of Occurrences');
xlim([0 xlimitMaxP]);
histogram(pmanFmag,'Normalization','probability','FaceColor','g','EdgeColor','g','BinWidth',binWdth);
plot(xFitP,pmanFmag_Y,'g','LineWidth',2);
histogram(p_ug_Fmag,'Normalization','probability','FaceColor','r','EdgeColor','r','FaceAlpha',0.01,'BinWidth',binWdth);
plot(xFitP,p_ug_Fmag_Y,'r','LineWidth',2);
histogram(p_g_Fmag,'Normalization','probability','FaceColor','b','EdgeColor','b','FaceAlpha',0.05,'BinWidth',binWdth);
legend('Manual','Unguided','Guided');
plot(xFitP,p_g_Fmag_Y,'b','LineWidth',2);

subplot(1,2,2); grid on; hold on;
xlabel('||Force|| (mN)'); ylabel('Percentage of Occurrences');
title('Cadaver Trials');
xlim([0 xlimitMaxC]);
histogram(cmanFmag,'Normalization','probability','FaceColor','g','EdgeColor','g','BinWidth',binWdth);
histogram(c_ug_Fmag,'Normalization','probability','FaceColor','r','EdgeColor','r','FaceAlpha',0.05,'BinWidth',binWdth);
histogram(c_g_Fmag,'Normalization','probability','FaceColor','b','EdgeColor','b','FaceAlpha',0.05,'BinWidth',binWdth);
legend('Manual','Unguided','Guided');