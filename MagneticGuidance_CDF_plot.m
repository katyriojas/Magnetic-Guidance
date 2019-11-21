% Katy Riojas
% Plot of CDF to compare all trials
% Last Updated: 11/21/19

% This plot needs to be using the trimmed data so that all (i.e., manual,
% phantom, and cadaver data are for a full insertion and that is it)

clc; addpath('functions');

regenerate_manual_data = false;
regenerate_phantom_data = false;
regenerate_cadaver_data = false;

if regenerate_manual_data
    LoadRALData_Manual; % regen data
elseif ~exist('data_manual_phantom','var')
    load('data\phantom\data_manual_phantom.mat'); % load already generated
elseif ~exist('data_manual_cadaver','var')
    load('data\cadaver\data_manual_cadaver.mat'); 
end

% guided/unguided phantom data
if regenerate_phantom_data
    LoadRALData_Robotic_Phantom % regen data
elseif ~exist('data_robotic_phantom','var') % if not already loaded
    load('data\phantom\data_robotic_phantom.mat'); % load already generated
end

if regenerate_cadaver_data
    LoadRALData_Robotic_Cadaver;
elseif ~exist('data_robotic_cadaver','var') % if not already loaded
    load('data\cadaver\data_robotic_cadaver.mat'); % load already generated
end

%% Plot variables
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

%% We first stack the data we care about into 3 vectors
phantom_m_Fmag = 0;
phantom_ug_Fmag = 0;
phantom_g_Fmag = 0;

for ii = 1:size(data_robotic_phantom,2)
    phantom_m_Fmag = [phantom_m_Fmag;data_manual_phantom(ii).Fmag_trimmed];
    phantom_ug_Fmag = [phantom_ug_Fmag;data_robotic_phantom(ii).nomag_mea.Fmag];
    phantom_g_Fmag = [phantom_g_Fmag;data_robotic_phantom(ii).mag.Fmag];
end

% Get rid of first element because that was our dummy zero
phantom_m_Fmag(1) = []; 
phantom_ug_Fmag(1) = []; 
phantom_g_Fmag(1) = [];

cadaver_m_Fmag = 0;
cadaver_ug_Fmag = 0;
cadaver_g_Fmag = 0;

for ii = 1:size(data_robotic_cadaver,2)
    cadaver_m_Fmag = [cadaver_m_Fmag;data_manual_cadaver(ii).Fmag_trimmed];
    cadaver_ug_Fmag = [cadaver_ug_Fmag;data_robotic_cadaver(ii).nomag_Fmagsmooth_trimmed];
    cadaver_g_Fmag = [cadaver_g_Fmag;data_robotic_cadaver(ii).mag_Fmagsmooth_trimmed];
end

% Get rid of first element because that was our dummy zero
cadaver_m_Fmag(1) = []; 
cadaver_ug_Fmag(1) = []; 
cadaver_g_Fmag(1) = [];

%% Plot CDF figure
figure(1); clf(1);
subplot(1,2,1); grid on; hold on;
xlabel('||Force|| (mN)');
title('(a) Phantom Trials');
ylabel({'Empirical CDF';'F(X)'});

[h1,X1,Y1,X1_inset,Y1_inset] = plotCDF(phantom_m_Fmag,colorMat(2,:),...
                                       LS,Fthreshold,upperThresh,lowerThresh);
[h2,X2,Y2,X2_inset,Y2_inset] = plotCDF(phantom_ug_Fmag,colorMat(1,:),...
                                       LS,Fthreshold,upperThresh,lowerThresh);
[h3,X3,Y3,X3_inset,Y3_inset] = plotCDF(phantom_g_Fmag,colorMat(3,:),...
                                       LS,Fthreshold,upperThresh,lowerThresh);

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
[h4,X4,Y4,X4_inset,Y4_inset] = plotCDF(cadaver_m_Fmag,colorMat(2,:),...
                                       LS,Fthreshold,upperThresh,lowerThresh);
[h5,X5,Y5,X5_inset,Y5_inset] = plotCDF(cadaver_ug_Fmag,colorMat(1,:),...
                                       LS,Fthreshold,upperThresh,lowerThresh);
[h6,X6,Y6,X6_inset,Y6_inset] = plotCDF(cadaver_g_Fmag,colorMat(3,:),...
                                       LS,Fthreshold,upperThresh,lowerThresh);

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