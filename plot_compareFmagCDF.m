% Katy Riojas
% 9/6/19
clear all; clc;
addpath('functions');

% Plot a histogram of Force Percentages during insertions
% Pull in manual data
load('data\phantom\pmanFmag.mat'); % load pmanFmag
load('data\cadaver\cmanFmag.mat'); % load cmanFmag

load('data\phantom\p_ug_Fmag.mat','p_ug_Fmag'); % load p_ug_Fmag
load('data\phantom\p_g_Fmag.mat','p_g_Fmag'); % load p_g_Fmag

load('data\cadaver\c_ug_Fmag.mat','c_ug_Fmag'); % load c_ug_Fmag
load('data\cadaver\c_g_Fmag.mat','c_g_Fmag'); % load c_g_Fmag

%% Plot CDF
colorMat = distinguishable_colors(3);
grey = [0.52,0.52,0.51];
buffer = 5;
buffer2 = 25;
LS = '-';
LS2 = '--';
Fthreshold = 25; %[mN]
upperThresh = 30; %[mN] upper and lower threshold used for inset
lowerThresh = 20; % [mN]
Msize = 30;
inset_width = 0.20;
inset_height = 0.34;
inset1_x = 0.22;
inset2_x = 0.68;
insety = 0.25;

figure(1); clf(1);
subplot(1,2,1); grid on; hold on;
xlabel('||Force|| (mN)');
title('(a) Phantom Trials');
ylabel({'Empirical CDF';'F(X)'});
[h1,X1,Y1,X1_inset,Y1_inset] = plotCDF(pmanFmag,colorMat(2,:),LS,Fthreshold,upperThresh,lowerThresh);
[h2,X2,Y2,X2_inset,Y2_inset] = plotCDF(p_ug_Fmag,colorMat(1,:),LS,Fthreshold,upperThresh,lowerThresh);
[h3,X3,Y3,X3_inset,Y3_inset] = plotCDF(p_g_Fmag,colorMat(3,:),LS,Fthreshold,upperThresh,lowerThresh);

% Draw horizontal lines:
line([0,X1],[Y1,Y1],'LineStyle',LS2,'LineWidth',1,'Color',colorMat(2,:));
line([0,X2],[Y2,Y2],'LineStyle',LS2,'LineWidth',1,'Color',colorMat(1,:));
line([0,X1],[Y3,Y3],'LineStyle',LS2,'LineWidth',1,'Color',colorMat(3,:));
% Draw vertical line:
line([Fthreshold,Fthreshold],[0,max([Y1,Y2,Y3])],...
    'LineStyle',LS2,'LineWidth',1,'Color',grey);
text(Fthreshold+buffer,0.05,'25mN');
text(Fthreshold+buffer2,Y1,strcat(num2str(Y1*100,'%2.0f'),'%'))
text(Fthreshold+buffer2,Y2,strcat(num2str(Y2*100,'%2.0f'),'%'))
text(Fthreshold+buffer2,Y3,strcat(num2str(Y3*100,'%2.0f'),'%'))
% Draw Markers:
scatter(X1,Y1,Msize,'filled','MarkerFaceColor',colorMat(2,:),'MarkerEdgeColor','None');
scatter(X2,Y2,Msize,'filled','MarkerFaceColor',colorMat(1,:),'MarkerEdgeColor','None');
scatter(X3,Y3,Msize,'filled','MarkerFaceColor',colorMat(3,:),'MarkerEdgeColor','None');
legend([h1,h2,h3],{'Case 1','Case 2','Case 3'});
% Draw Inset
axes('position',[inset1_x,insety,inset_width,inset_height]);
grid on; hold on; ylim([0.4 0.9]);
box on; axis tight; % put box around new pair of axes
xlabel('||Force|| (mN)');
plot(X1_inset,Y1_inset,'Color',colorMat(2,:),'LineStyle','-','LineWidth',1);
plot(X2_inset,Y2_inset,'Color',colorMat(1,:),'LineStyle','-','LineWidth',1);
plot(X3_inset,Y3_inset,'Color',colorMat(3,:),'LineStyle','-','LineWidth',1);
scatter(X1,Y1,Msize,'filled','MarkerFaceColor',colorMat(2,:),'MarkerEdgeColor','None');
scatter(X2,Y2,Msize,'filled','MarkerFaceColor',colorMat(1,:),'MarkerEdgeColor','None');
scatter(X3,Y3,Msize,'filled','MarkerFaceColor',colorMat(3,:),'MarkerEdgeColor','None');

% Cadaver Data
subplot(1,2,2); grid on; hold on;
title('(b) Cadaver Trials');
xlabel('||Force|| (mN)');
[h4,X4,Y4,X4_inset,Y4_inset] = plotCDF(cmanFmag,colorMat(2,:),LS,Fthreshold,upperThresh,lowerThresh);
[h5,X5,Y5,X5_inset,Y5_inset] = plotCDF(c_ug_Fmag,colorMat(1,:),LS,Fthreshold,upperThresh,lowerThresh);
[h6,X6,Y6,X6_inset,Y6_inset] = plotCDF(c_g_Fmag,colorMat(3,:),LS,Fthreshold,upperThresh,lowerThresh);

% Draw horizontal lines:
line([0,X4],[Y4,Y4],'LineStyle',LS2,'LineWidth',1,'Color',colorMat(2,:));
line([0,X5],[Y5,Y5],'LineStyle',LS2,'LineWidth',1,'Color',colorMat(1,:));
line([0,X6],[Y6,Y6],'LineStyle',LS2,'LineWidth',1,'Color',colorMat(3,:));
% Draw vertical line:
line([Fthreshold,Fthreshold],[0,max([Y4,Y5,Y6])],...
    'LineStyle',LS2,'LineWidth',1,'Color',grey);
text(Fthreshold+buffer,0.05,'25mN')
text(Fthreshold+buffer2,Y4,strcat(num2str(Y4*100,'%2.0f'),'%'))
text(Fthreshold+buffer2,Y5,strcat(num2str(Y5*100,'%2.0f'),'%'))
text(Fthreshold+buffer2,Y6,strcat(num2str(Y6*100,'%2.0f'),'%'))
% Draw Markers
scatter(X4,Y4,Msize,'filled','MarkerFaceColor',colorMat(2,:),'MarkerEdgeColor','None');
scatter(X5,Y5,Msize,'filled','MarkerFaceColor',colorMat(1,:),'MarkerEdgeColor','None');
scatter(X6,Y6,Msize,'filled','MarkerFaceColor',colorMat(3,:),'MarkerEdgeColor','None');
%Legend
legend([h4,h5,h6],{'Case 1','Case 2','Case 3'});

% Draw Inset
axes('position',[inset2_x,insety,inset_width,inset_height]);
grid on; hold on; ylim([0.4 0.9]);
box on; axis tight; % put box around new pair of axes
xlabel('||Force|| (mN)');
plot(X4_inset,Y4_inset,'Color',colorMat(2,:),'LineStyle','-','LineWidth',1);
plot(X5_inset,Y5_inset,'Color',colorMat(1,:),'LineStyle','-','LineWidth',1);
plot(X6_inset,Y6_inset,'Color',colorMat(3,:),'LineStyle','-','LineWidth',1);
scatter(X4,Y4,Msize,'filled','MarkerFaceColor',colorMat(2,:),'MarkerEdgeColor','None');
scatter(X5,Y5,Msize,'filled','MarkerFaceColor',colorMat(1,:),'MarkerEdgeColor','None');
scatter(X6,Y6,Msize,'filled','MarkerFaceColor',colorMat(3,:),'MarkerEdgeColor','None');


matPerc = [Y1,Y2,Y3;...
            Y4,Y5,Y6]
% %% Plot PDF
% % numbins = 50;
% binWdth = 1;
% xlimitMaxP = max([pmanFmag;p_ug_Fmag;p_g_Fmag]);
% xFitP = linspace(0,xlimitMaxP,500);
% xlimitMaxC = max([cmanFmag;c_ug_Fmag;c_g_Fmag]);
% xFitC = linspace(0,xlimitMaxC,500);
% 
% % Fit gamma
% pmanFmag_dist = fitdist(sort(pmanFmag),'gamma');
% pmanFmag_dist.ParameterValues
% pmanFmag_Y = gampdf(xFitP,pmanFmag_dist.ParameterValues(1),pmanFmag_dist.ParameterValues(2));
% 
% p_ug_Fmag_dist = fitdist(sort(p_ug_Fmag),'gamma');
% p_ug_Fmag_dist.ParameterValues
% p_ug_Fmag_Y = gampdf(xFitP,p_ug_Fmag_dist.ParameterValues(1),p_ug_Fmag_dist.ParameterValues(2));
% 
% p_g_Fmag_dist = fitdist(sort(p_g_Fmag),'gamma');
% p_g_Fmag_dist.ParameterValues
% p_g_Fmag_Y = gampdf(xFitP,p_g_Fmag_dist.ParameterValues(1),p_g_Fmag_dist.ParameterValues(2));
% 
% figure(2); clf(2);
% subplot(1,2,1); grid on; hold on;
% title('Phantom Insertions');
% xlabel('||Force|| (mN)'); ylabel('Percentage of Occurrences');
% xlim([0 xlimitMaxP]);
% histogram(pmanFmag,'Normalization','probability','FaceColor','g','EdgeColor','g','BinWidth',binWdth);
% plot(xFitP,pmanFmag_Y,'g','LineWidth',2);
% histogram(p_ug_Fmag,'Normalization','probability','FaceColor','r','EdgeColor','r','FaceAlpha',0.01,'BinWidth',binWdth);
% plot(xFitP,p_ug_Fmag_Y,'r','LineWidth',2);
% histogram(p_g_Fmag,'Normalization','probability','FaceColor','b','EdgeColor','b','FaceAlpha',0.05,'BinWidth',binWdth);
% legend('Manual','Unguided','Guided');
% plot(xFitP,p_g_Fmag_Y,'b','LineWidth',2);
% 
% subplot(1,2,2); grid on; hold on;
% xlabel('||Force|| (mN)'); ylabel('Percentage of Occurrences');
% title('Cadaver Trials');
% xlim([0 xlimitMaxC]);
% histogram(cmanFmag,'Normalization','probability','FaceColor','g','EdgeColor','g','BinWidth',binWdth);
% histogram(c_ug_Fmag,'Normalization','probability','FaceColor','r','EdgeColor','r','FaceAlpha',0.05,'BinWidth',binWdth);
% histogram(c_g_Fmag,'Normalization','probability','FaceColor','b','EdgeColor','b','FaceAlpha',0.05,'BinWidth',binWdth);
% legend('Manual','Unguided','Guided');