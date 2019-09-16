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

%% Pull in phantom AID data
xlabel('Insertion Depth (deg)');
ylabel('||Force|| (mN)');
load('data\phantom\manual\pman1_angular_depth.mat');
pman1_angular_depth = insertion_angle;
load('data\phantom\manual\pman2_angular_depth.mat');
pman2_angular_depth = insertion_angle;
load('data\phantom\manual\pman3_angular_depth.mat');
pman3_angular_depth = insertion_angle;
load('data\phantom\manual\pman4_angular_depth.mat');
pman4_angular_depth = insertion_angle;

% Interpolate force points
% pman1.interp_angdepth = interp1(pman1_angular_depth.time,...
%     pman1_angular_depth.angle, data_pman1.time(p1_vec)-...
%     data_pman1.time(p1_vec(1)), 'linear', 'extrap');
% pman2.interp_angdepth = interp1(pman2_angular_depth.time,...
%     pman2_angular_depth.angle, data_pman2.time(p2_vec)-...
%     data_pman2.time(p2_vec(1)), 'linear', 'extrap');
% pman3.interp_angdepth = interp1(pman3_angular_depth.time,...
%     pman3_angular_depth.angle, data_pman3.time(p3_vec)-...
%     data_pman3.time(p3_vec(1)), 'linear', 'extrap');
% pman4.interp_angdepth = interp1(pman4_angular_depth.time,...
%     pman4_angular_depth.angle, data_pman4.time(p4_vec)-...
%     data_pman4.time(p4_vec(1)), 'linear', 'extrap');
pman1.interp_angdepth = interp1(pman1_angular_depth.time,...
    pman1_angular_depth.angle, data_pman1.time(p1_vec)-...
    data_pman1.time(p1_vec(1)), 'linear', 'extrap');
pman2.interp_angdepth = interp1(pman2_angular_depth.time,...
    pman2_angular_depth.angle, data_pman2.time(p2_vec)-...
    data_pman2.time(p2_vec(1)), 'linear', 'extrap');
pman3.interp_angdepth = interp1(pman3_angular_depth.time,...
    pman3_angular_depth.angle, data_pman3.time(p3_vec)-...
    data_pman3.time(p3_vec(1)), 'linear', 'extrap');
pman4.interp_angdepth = interp1(pman4_angular_depth.time,...
    pman4_angular_depth.angle, data_pman4.time(p4_vec)-...
    data_pman4.time(p4_vec(1)), 'linear', 'extrap');

%% Plot 1: Averages
% if exist('data','var')
load('data\phantom\data From Trevor\data_phantom.mat');
% end
%%
% Crop Fmag_smooth to the start and end depths
Fmag_smooth_man1 = data_pman1.Fmag_smooth(p1_vec);
Fmag_smooth_man2 = data_pman2.Fmag_smooth(p2_vec);
Fmag_smooth_man3 = data_pman3.Fmag_smooth(p3_vec);
Fmag_smooth_man4 = data_pman4.Fmag_smooth(p4_vec);

pman1_angdepth = pman1.interp_angdepth;
pman2_angdepth = pman2.interp_angdepth;
pman3_angdepth = pman3.interp_angdepth;
pman4_angdepth = pman4.interp_angdepth;
% save('ManualPhantomData.mat','pman1_angdepth','pman2_angdepth',...
%     'pman3_angdepth','pman4_angdepth','Fmag_smooth_man1',...
%     'Fmag_smooth_man2','Fmag_smooth_man3','Fmag_smooth_man4');

Fmag_smooth_ug1 = data(1).nomag_mea.Fmag_smooth;
Fmag_smooth_ug2 = data(2).nomag_mea.Fmag_smooth;
Fmag_smooth_ug3 = data(3).nomag_mea.Fmag_smooth;
Fmag_smooth_ug4 = data(4).nomag_mea.Fmag_smooth;

Fmag_smooth_g1 = data(1).mag.Fmag_smooth;
Fmag_smooth_g2 = data(2).mag.Fmag_smooth;
Fmag_smooth_g3 = data(3).mag.Fmag_smooth;
Fmag_smooth_g4 = data(4).mag.Fmag_smooth;

% Variables for averaging
degspan = 1;
maxX = ceil(max([pman1.interp_angdepth;...
                 pman2.interp_angdepth;...
                 pman3.interp_angdepth;...
                 pman4.interp_angdepth]));
             
numsteps = ceil(maxX/degspan);           
stepsize = 1; 
AOIvec = linspace(0,maxX,numsteps);

Fmag_avg_man = zeros(length(AOIvec),1);
pman_std = zeros(length(AOIvec),1);

Fmag_avg_ug = zeros(size(AOIvec));
ug_std = zeros(size(AOIvec));

Fmag_avg_g = zeros(size(AOIvec));
g_std = zeros(size(AOIvec));

Fmag_man1_bin = zeros(size(AOIvec));
Fmag_man2_bin = zeros(size(AOIvec));
Fmag_man3_bin = zeros(size(AOIvec));
Fmag_man4_bin = zeros(size(AOIvec));

for ii = 2:length(AOIvec)
    
    idx1 = find((pman1.interp_angdepth<AOIvec(ii))&(pman1.interp_angdepth>AOIvec(ii-1)));
    idx2 = find((pman2.interp_angdepth<AOIvec(ii))&(pman2.interp_angdepth>AOIvec(ii-1)));
    idx3 = find((pman3.interp_angdepth<AOIvec(ii))&(pman3.interp_angdepth>AOIvec(ii-1)));
    idx4 = find((pman4.interp_angdepth<AOIvec(ii))&(pman4.interp_angdepth>AOIvec(ii-1)));
    
    Fmag_man1 = Fmag_smooth_man1(idx1);
    Fmag_man2 = Fmag_smooth_man2(idx2);
    Fmag_man3 = Fmag_smooth_man3(idx3);
    Fmag_man4 = Fmag_smooth_man3(idx4);
    
    Fmag_man1_bin(ii) = mean(Fmag_man1);
    Fmag_man2_bin(ii) = mean(Fmag_man2);
    Fmag_man3_bin(ii) = mean(Fmag_man3);
    Fmag_man4_bin(ii) = mean(Fmag_man4);

end

man1_binned = [AOIvec',Fmag_man1_bin'];
man1_binned(any(isnan(man1_binned),2),:) = [];

man2_binned = [AOIvec',Fmag_man2_bin'];
man2_binned(any(isnan(man2_binned),2),:) = [];

man3_binned = [AOIvec',Fmag_man3_bin'];
man3_binned(any(isnan(man3_binned),2),:) = [];

man4_binned = [AOIvec',Fmag_man4_bin'];
man4_binned(any(isnan(man4_binned),2),:) = [];

pman(1).angdepth = man1_binned(:,1);
pman(1).Fmag_smooth = man1_binned(:,2);

pman(2).angdepth = man2_binned(:,1);
pman(2).Fmag_smooth = man2_binned(:,2);

pman(3).angdepth = man3_binned(:,1);
pman(3).Fmag_smooth = man3_binned(:,2);

pman(4).angdepth = man4_binned(:,1);
pman(4).Fmag_smooth = man4_binned(:,2);

% rough plot
figure; clf; grid on; hold on;
for ii = 1:4
    plot(pman(ii).angdepth,pman(ii).Fmag_smooth);
end
save('ManualPhantomData_updated.mat','pman');

for ii = 2:length(AOIvec)
    
    idx1 = find((pman1.interp_angdepth<AOIvec(ii))&(pman1.interp_angdepth>AOIvec(ii-1)));
    idx2 = find((pman2.interp_angdepth<AOIvec(ii))&(pman2.interp_angdepth>AOIvec(ii-1)));
    idx3 = find((pman3.interp_angdepth<AOIvec(ii))&(pman3.interp_angdepth>AOIvec(ii-1)));
    idx4 = find((pman4.interp_angdepth<AOIvec(ii))&(pman4.interp_angdepth>AOIvec(ii-1)));
    
    Fmag_smooth_bin = [Fmag_smooth_man1(idx1);...
                       Fmag_smooth_man2(idx2);...
                       Fmag_smooth_man3(idx3);...
                       Fmag_smooth_man4(idx4)];
                
    Fmag_avg_man(ii) = mean(Fmag_smooth_bin);
    pman_std(ii) = std(Fmag_smooth_bin);
    
    idx1 = find((data(1).nomag_mea_interp_angdepth<AOIvec(ii))...
                &(data(1).nomag_mea_interp_angdepth>AOIvec(ii-1)));
    idx2 = find((data(2).nomag_mea_interp_angdepth<AOIvec(ii))...
                &(data(2).nomag_mea_interp_angdepth>AOIvec(ii-1)));
    idx3 = find((data(3).mag_interp_angdepth<AOIvec(ii))...
                &(data(3).mag_interp_angdepth>AOIvec(ii-1)));
    idx4 = find((data(4).nomag_mea_interp_angdepth<AOIvec(ii))...
                &(data(4).nomag_mea_interp_angdepth>AOIvec(ii-1)));
    
    Fmag_smooth_bin = [Fmag_smooth_ug1(idx1);...
                       Fmag_smooth_ug2(idx2);...
                       Fmag_smooth_man3(idx3);
                       Fmag_smooth_ug4(idx4)];
                   
    Fmag_avg_ug(ii) = mean(Fmag_smooth_bin);
    ug_std(ii) = std(Fmag_smooth_bin);
    
    idx1 = find((data(1).mag_interp_angdepth<AOIvec(ii))...
                &(data(1).mag_interp_angdepth>AOIvec(ii-1)));
    idx2 = find((data(2).mag_interp_angdepth<AOIvec(ii))...
                &(data(2).mag_interp_angdepth>AOIvec(ii-1)));
    idx3 = find((data(3).mag_interp_angdepth<AOIvec(ii))...
                &(data(3).mag_interp_angdepth>AOIvec(ii-1)));
    idx4 = find((data(4).mag_interp_angdepth<AOIvec(ii))...
                &(data(4).mag_interp_angdepth>AOIvec(ii-1)));
    
    Fmag_smooth_bin = [Fmag_smooth_g1(idx1);...
                       Fmag_smooth_g2(idx2);...
                       Fmag_smooth_g3(idx3);
                       Fmag_smooth_g4(idx4)];
                   
    Fmag_avg_g(ii) = mean(Fmag_smooth_bin);
    g_std(ii) = std(Fmag_smooth_bin);

end

% Check manual data fit
figure(1); clf(1); grid on; hold on;
xlim([0,maxX]); 
title('Case 1 Data');
xlabel('Angular Insertion Depth [deg]'); ylabel('Average ||Force|| [mN]');
plot(pman1.interp_angdepth,data_pman1.Fmag_smooth(p1_vec),'Color',...
    cMat(1,:),'LineWidth',1,'LineStyle',':');
plot(pman2.interp_angdepth,data_pman2.Fmag_smooth(p2_vec),'Color',...
    cMat(2,:),'LineWidth',1,'LineStyle',':');
plot(pman3.interp_angdepth,data_pman3.Fmag_smooth(p3_vec),'Color',...
    cMat(3,:),'LineWidth',1,'LineStyle',':');
plot(pman4.interp_angdepth,data_pman4.Fmag_smooth(p4_vec),'Color',...
    cMat(4,:),'LineWidth',1,'LineStyle',':');
plot(AOIvec-degspan/2, Fmag_avg_man, 'Color', 'k','LineWidth',1,'LineStyle','-');

% Check no mag data fit
figure(2); clf(2); grid on; hold on;
title('Case 2 Data');
xlabel('Angular Insertion Depth [deg]'); ylabel('Average ||Force|| [mN]');
plot(data(1).nomag_mea_interp_angdepth,Fmag_smooth_ug1,'Color',...
    cMat(1,:),'LineWidth',1,'LineStyle',':');
plot(data(2).nomag_mea_interp_angdepth,Fmag_smooth_ug2,'Color',...
    cMat(2,:),'LineWidth',1,'LineStyle',':');
plot(data(3).nomag_mea_interp_angdepth,Fmag_smooth_ug3,'Color',...
    cMat(3,:),'LineWidth',1,'LineStyle',':');
plot(data(4).nomag_mea_interp_angdepth,Fmag_smooth_ug4,'Color',...
    cMat(4,:),'LineWidth',1,'LineStyle',':');
plot(AOIvec-degspan/2, Fmag_avg_ug,...
    'Color', 'k','LineWidth',1,'LineStyle','-');

% Check mag data fit
figure(3); clf(3); grid on; hold on;
title('Case 3 Data');
xlabel('Angular Insertion Depth [deg]'); ylabel('Average ||Force|| [mN]');
plot(data(1).mag_interp_angdepth,Fmag_smooth_g1,'Color',...
    cMat(1,:),'LineWidth',1,'LineStyle',':');
plot(data(2).mag_interp_angdepth,Fmag_smooth_g2,'Color',...
    cMat(2,:),'LineWidth',1,'LineStyle',':');
plot(data(3).mag_interp_angdepth,Fmag_smooth_g3,'Color',...
    cMat(3,:),'LineWidth',1,'LineStyle',':');
plot(data(4).mag_interp_angdepth,Fmag_smooth_g4,'Color',...
    cMat(4,:),'LineWidth',1,'LineStyle',':');
plot(AOIvec-degspan/2, Fmag_avg_g,...
    'Color', 'k','LineWidth',1,'LineStyle','-');

% Now Plot them all!
figure(4); clf(4); grid on; hold on; xlim([0 400]);
xlabel('Angular Insertion Depth (\circ)'); ylabel('Average ||Force|| (mN)');
plot(AOIvec-degspan/2, Fmag_avg_man, 'Color', 'r','LineWidth',1,'LineStyle','-');
plot(AOIvec-degspan/2, Fmag_avg_ug, 'Color', 'b','LineWidth',1,'LineStyle','-');
plot(AOIvec-degspan/2, Fmag_avg_g, 'Color', 'g','LineWidth',1,'LineStyle','-');
legend('Case 1','Case 2','Case 3');
% fill([(AOIvec-degspan/2)' fliplr((AOIvec-degspan/2)')],[Fmag_smooth_avg + Fmag_smooth_std,...
%     fliplr(Fmag_smooth_avg - Fmag_smooth_std)],'k','FaceAlpha',0.2,...
%     'LineStyle','none');