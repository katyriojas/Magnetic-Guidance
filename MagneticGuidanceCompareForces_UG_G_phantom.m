%% Import force/smaract data & force sensor calibration

% clear all; close all; clc;

addpath('classDef','functions','data');
load('avg_cal_slopes.mat'); % loads cal_slopes


% Filepaths to CSVs exported from ROS bags

base_path = 'data\phantom\';

% unguided, unmodified EA
filepaths(1).nomag_ea.smaract = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial1_1.25_smaract.csv');
filepaths(1).nomag_ea.force   = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial1_1.25_force.csv');

filepaths(2).nomag_ea.smaract = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial2_1.25_smaract.csv');
filepaths(2).nomag_ea.force   = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial2_1.25_force.csv');

filepaths(3).nomag_ea.smaract = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial3_1.25_smaract.csv');
filepaths(3).nomag_ea.force   = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial3_1.25_force.csv');

filepaths(4).nomag_ea.smaract = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial4_1.25_smaract.csv');
filepaths(4).nomag_ea.force   = fullfile(base_path, 'unguided_ea_saline_1.25\phantom_ug_ea2_trial4_1.25_force.csv');

% unguided, magnet-tipped EA
filepaths(1).nomag_mea.smaract = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial1_1.25_smaract.csv');
filepaths(1).nomag_mea.force   = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial1_1.25_force.csv');

filepaths(2).nomag_mea.smaract = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial2_1.25_smaract.csv');
filepaths(2).nomag_mea.force   = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial2_1.25_force.csv');

filepaths(3).nomag_mea.smaract = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial3_1.25_smaract.csv');
filepaths(3).nomag_mea.force   = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial3_1.25_force.csv');

filepaths(4).nomag_mea.smaract = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial4_1.25_smaract.csv');
filepaths(4).nomag_mea.force   = fullfile(base_path, 'unguided_mea_saline_1.25\phantom_ug_mea1_trial4_1.25_force.csv');

% guided (with magnet)
filepaths(1).mag.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial1_1.25_smaract.csv');
filepaths(1).mag.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial1_1.25_force.csv');

filepaths(2).mag.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial2_1.25_smaract.csv');
filepaths(2).mag.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial2_1.25_force.csv');

filepaths(3).mag.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial3_1.25_smaract.csv');
filepaths(3).mag.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial3_1.25_force.csv');

filepaths(4).mag.smaract = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial4_1.25_smaract.csv');
filepaths(4).mag.force   = fullfile(base_path, 'guided_saline_1.25\phantom_g_mea1_trial4_1.25_force.csv');


%% Create data objects

for ii=1:length(filepaths)
    data(ii).nomag_ea  = MagneticGuidanceData(filepaths(ii).nomag_ea,  cal_slopes);
    data(ii).nomag_mea = MagneticGuidanceData(filepaths(ii).nomag_mea, cal_slopes);
    data(ii).mag       = MagneticGuidanceData(filepaths(ii).mag,       cal_slopes);
end


%% set smoothing span (proportion of data points)
force_smooth_span = 40; % [# samples]

for ii=1:length(filepaths)
    data(ii).nomag_ea.smooth_span  = force_smooth_span;
    data(ii).nomag_mea.smooth_span = force_smooth_span;
    data(ii).mag.smooth_span       = force_smooth_span;
end


%% determine angular insertion depths from processed video

angle_smooth_span = 20;

base_path = 'D:\Trevor\My Documents\MED lab\Cochlear R01\Mag Steering\Experiments\RAL\';

filepaths(1).nomag_mea.vid = fullfile(base_path, 'phantom_ug_mea1_trial1_1.25\trial1-ug-mea1-tracked.mp4');
filepaths(1).mag.vid       = fullfile(base_path, 'phantom_g_mea1_trial1_1.25\trial1-guided-mea1-1.25-tracked.mp4');

filepaths(2).nomag_mea.vid = fullfile(base_path, 'phantom_ug_mea1_trial2_1.25\trial2-ug-mea1-tracked.mp4');
filepaths(2).mag.vid       = fullfile(base_path, 'phantom_g_mea1_trial2_1.25\trial2-g-mea1-tracked.mp4');

filepaths(3).nomag_mea.vid = fullfile(base_path, 'phantom_ug_mea1_trial3_1.25\trial3-ug-mea1-tracked.mp4');
filepaths(3).mag.vid       = fullfile(base_path, 'phantom_g_mea1_trial3_1.25\trial3-g-mea1-tracked.mp4');

filepaths(4).nomag_mea.vid = fullfile(base_path, 'phantom_ug_mea1_trial4_1.25\trial4-ug-mea1-tracked.mp4');
filepaths(4).mag.vid       = fullfile(base_path, 'phantom_g_mea1_trial4_1.25\trial4-g-mea1-tracked.mp4');

% for ii = 1:2
for ii = 1:length(filepaths)
    
    % segment video frames to determine angular depth
    if isempty(data(ii).mag_angular_depth) % only compute if needed
        data(ii).mag_angular_depth = MagneticGuidanceGetAngleFromVideo(filepaths(ii).mag.vid, angle_smooth_span);
    end

    if isempty(data(ii).nomag_mea_angular_depth)
        data(ii).nomag_mea_angular_depth = MagneticGuidanceGetAngleFromVideo(filepaths(ii).nomag_mea.vid, angle_smooth_span);
    end
    


    % interpolate to find angular depth at each force measurement time
    if isempty(data(ii).mag_interp_angdepth)
        data(ii).mag_interp_angdepth       = interp1(data(ii).mag_angular_depth.time,...
                                                     data(ii).mag_angular_depth.angle_smooth,...
                                                     data(ii).mag.time_insertion, 'linear', 'extrap');
    end

    if isempty(data(ii).nomag_mea_interp_angdepth) && ii~=3 % video stopped early on nomag trial 3  
        data(ii).nomag_mea_interp_angdepth = interp1(data(ii).nomag_mea_angular_depth.time,...
                                                     data(ii).nomag_mea_angular_depth.angle_smooth,...
                                                     data(ii).nomag_mea.time_insertion, 'linear', 'extrap');
    end

end

%% use nomag trial 4 to fill remainder of nomag trial 3

% trial 3 data until stop point
nomag3_vidstop_time  = data(3).nomag_mea_angular_depth.time(end);
nomag3_stop_ind = find(data(3).nomag_mea.time_insertion > nomag3_vidstop_time, 1);

data(3).nomag_mea_interp_angdepth = interp1(data(3).nomag_mea_angular_depth.time,...
                                            data(3).nomag_mea_angular_depth.angle_smooth,...
                                            data(3).nomag_mea.time_insertion(1:nomag3_stop_ind), 'linear', 'extrap');

% addition from trial 4
nomag3_vidstop_depth = data(3).nomag_mea.depth_insertion( find(data(3).nomag_mea.time_insertion > nomag3_vidstop_time, 1) );
nomag3_total_time = data(3).nomag_mea.time_insertion(end);
nomag3_total_depth = data(3).nomag_mea.depth_insertion(end);


nomag4_inds = find(data(4).nomag_mea.depth_insertion > nomag3_vidstop_depth, 1): find(data(4).nomag_mea.depth_insertion > nomag3_total_depth, 1);

nomag3_interp_angdepth = [data(3).nomag_mea_interp_angdepth; data(4).nomag_mea_interp_angdepth(nomag4_inds)];

nomag3_insertion_depths = linspace(data(3).nomag_mea.depth_insertion(1), data(3).nomag_mea.depth_insertion(end), length(data(3).nomag_mea_interp_angdepth));

% re-interpolate
data(3).nomag_mea_interp_angdepth = interp1(nomag3_insertion_depths,...
                                            data(3).nomag_mea_interp_angdepth,...
                                            data(3).nomag_mea.depth_insertion, 'linear', 'extrap');



%% Plot -> linear vs angular insertion depth

figure(25); clf(25); 
hold on; grid on;

colors = distinguishable_colors(2*length(data)+1);

% for ii=3:4
for ii=1:length(data)
    h_mag(ii)   = plot( data(ii).mag.depth_insertion, data(ii).mag_interp_angdepth, 'Color',colors(ii*2-1,:), 'LineWidth', 2);
    h_nomag(ii) = plot( data(ii).nomag_mea.depth_insertion, data(ii).nomag_mea_interp_angdepth, 'Color',colors(ii*2,:), 'LineStyle',':', 'LineWidth', 2);
end

xlabel('Actuator Insertion Distance (mm)')
ylabel('Angular Insertion Depth (deg)')

labels.mag = {'mag1','mag2','mag3','mag4'};
labels.nomag = {'nomag1','nomag2','nomag3','nomag4'};
legend([h_mag(1:ii), h_nomag(1:ii)], [labels.mag(1:ii), labels.nomag(1:ii)], 'Location','nw')


%% Plot -> angular insertion depth vs force

% force components
figure(1); clf(1); 
hold on; grid on;

% for ii=1:2
for ii=1:length(data)

                   plot( data(ii).mag_interp_angdepth, data(ii).mag.Fx,        'Color', [1,0,0, 0.3], 'LineWidth', 1);
    h(ii).mag(1) = plot( data(ii).mag_interp_angdepth, data(ii).mag.Fx_smooth, 'Color', [1,0,0, 1.0], 'LineWidth', 2);
                   plot( data(ii).mag_interp_angdepth, data(ii).mag.Fy,        'Color', [0,1,0, 0.3], 'LineWidth', 1);
    h(ii).mag(2) = plot( data(ii).mag_interp_angdepth, data(ii).mag.Fy_smooth, 'Color', [0,1,0, 1.0], 'LineWidth', 2);
                   plot( data(ii).mag_interp_angdepth, data(ii).mag.Fz,        'Color', [0,0,1, 0.3], 'LineWidth', 1);
    h(ii).mag(3) = plot( data(ii).mag_interp_angdepth, data(ii).mag.Fz_smooth, 'Color', [0,0,1, 1.0], 'LineWidth', 2);

                     plot( data(ii).nomag_mea_interp_angdepth, data(ii).nomag_mea.Fx,        'Color', [0.95,0.21,0.62, 0.3], 'LineWidth', 1);
    h(ii).nomag(1) = plot( data(ii).nomag_mea_interp_angdepth, data(ii).nomag_mea.Fx_smooth, 'Color', [0.95,0.21,0.62, 1.0], 'LineWidth', 2);
                     plot( data(ii).nomag_mea_interp_angdepth, data(ii).nomag_mea.Fy,        'Color', [0.75,0.88,0.24, 0.3], 'LineWidth', 1);
    h(ii).nomag(2) = plot( data(ii).nomag_mea_interp_angdepth, data(ii).nomag_mea.Fy_smooth, 'Color', [0.75,0.88,0.24, 1.0], 'LineWidth', 2);
                     plot( data(ii).nomag_mea_interp_angdepth, data(ii).nomag_mea.Fz,        'Color', [0.21,0.65,0.95, 0.3], 'LineWidth', 1);
    h(ii).nomag(3) = plot( data(ii).nomag_mea_interp_angdepth, data(ii).nomag_mea.Fz_smooth, 'Color', [0.21,0.65,0.95, 1.0], 'LineWidth', 2);
end

title('Force vs Angular Insertion Depth')
xlabel('angular insertion depth (deg)')
ylabel('force (mN)')
legend(h(1).mag, {'Fx', 'Fy', 'Fz'}, 'Location', 'sw')


% force magnitude
figure(2); clf(2);
hold on; grid on;

for ii=1:2
% for ii=1:length(data)
    plot( data(ii).mag_interp_angdepth, data(ii).mag.Fmag,        'Color', [colors(ii*2-1,:), 0.3], 'LineWidth', 1);
    plot( data(ii).mag_interp_angdepth, data(ii).mag.Fmag_smooth, 'Color', [colors(ii*2-1,:), 1], 'LineWidth', 2);
    plot( data(ii).nomag_mea_interp_angdepth, data(ii).nomag_mea.Fmag,        'Color', [colors(ii*2,:), 0.3], 'LineWidth', 1, 'LineStyle',':', 'LineWidth', 1);
    plot( data(ii).nomag_mea_interp_angdepth, data(ii).nomag_mea.Fmag_smooth, 'Color', [colors(ii*2,:), 1],   'LineWidth', 1, 'LineStyle',':', 'LineWidth', 2);
end

title('Force vs Angular Insertion Depth')
xlabel('angular insertion depth (deg)')
ylabel('force (mN)')


%% plot -> linear insertion depth vs Fmag

alpha = 1; % reduce transparency of unguided plot lines

colorsMat = distinguishable_colors(3*length(data));

% Figure: force magnitude
figure(4); clf(4); hold on; grid on;


for ii=1:length(data)

                    plot( data(ii).nomag_ea.depth_insertion, data(ii).nomag_ea.Fmag,    'Color',[colorsMat(ii,:), 0.3*alpha],      'LineWidth',1);
  h(ii).nomag_ea =  plot( data(ii).nomag_ea.depth_insertion, data(ii).nomag_ea.Fmag_smooth,'Color', colorsMat(ii,:), 'LineStyle',':', 'LineWidth',1.5);
 
                    plot( data(ii).nomag_mea.depth_insertion, data(ii).nomag_mea.Fmag,        'Color', [colorsMat(ii*2,:),  0.3*alpha],      'LineWidth',1);
  h(ii).nomag_mea = plot( data(ii).nomag_mea.depth_insertion, data(ii).nomag_mea.Fmag_smooth, 'Color',  colorsMat(ii*2,:), 'LineStyle','--', 'LineWidth',1.5);

                    plot( data(ii).mag.depth_insertion, data(ii).mag.Fmag,       'Color',[colorsMat(ii*3,:), 0.3*alpha],  'LineWidth',1);
  h(ii).mag =       plot( data(ii).mag.depth_insertion, data(ii).mag.Fmag_smooth,'Color', colorsMat(ii*3,:), 'LineWidth',1.5);

end

title('Guided vs Unguided Insertion')
xlabel('Insertion Depth [mm]')
ylabel('||Force|| [mN]')

legend([h.nomag_ea, h.nomag_mea, h.mag],{'nomag1-ea', 'nomag2-ea','nomag3-ea','nomag4-ea',...
                                         'nomag1-mea','nomag2-mea','nomag3-mea','nomag4-mea',...
                                         'mag1','mag2','mag3','mag4'});


%% TO DO

xvec = linspace(0.25,23.7,1000);
% xvec = linspace(0.25,26.4,1000);

Fnomag_mea = [interp1(data_nomag1_mea.depth_insertion(2:end-1),data_nomag1_mea.Fmag_smooth(2:end-1),xvec);...
              interp1(data_nomag2_mea.depth_insertion(2:end-1),data_nomag2_mea.Fmag_smooth(2:end-1),xvec);...
              interp1(data_nomag3_mea.depth_insertion(2:end-1),data_nomag3_mea.Fmag_smooth(2:end-1),xvec);...
              interp1(data_nomag4_mea.depth_insertion(2:end-1),data_nomag4_mea.Fmag_smooth(2:end-1),xvec)];
% 
Fmag = [interp1(data_mag1.depth_insertion(2:end-1),data_mag1.Fmag_smooth(2:end-1),xvec);...
        interp1(data_mag2.depth_insertion(2:end-1),data_mag2.Fmag_smooth(2:end-1),xvec);...
        interp1(data_mag3.depth_insertion(2:end-1),data_mag3.Fmag_smooth(2:end-1),xvec);...
        interp1(data_mag4.depth_insertion(2:end-1),data_mag4.Fmag_smooth(2:end-1),xvec)];
    
Fnomag_ea = [interp1(data_nomag1_ea.depth_insertion(2:end-1),data_nomag1_ea.Fmag_smooth(2:end-1),xvec);...
             interp1(data_nomag2_ea.depth_insertion(2:end-1),data_nomag2_ea.Fmag_smooth(2:end-1),xvec);...
             interp1(data_nomag3_ea.depth_insertion(2:end-1),data_nomag3_ea.Fmag_smooth(2:end-1),xvec);...
             interp1(data_nomag4_ea.depth_insertion(2:end-1),data_nomag4_ea.Fmag_smooth(2:end-1),xvec)];

Favg_nomag_mea = mean(Fnomag_mea,1);
std_nomag_mea = std(Fnomag_mea);

Favg_mag = mean(Fmag,1);
std_mag = std(Fmag);

Favg_nomag_ea = nanmean(Fnomag_ea,1);
std_nomag_ea = std(Fnomag_ea);

colorsMat2 = distinguishable_colors(6);

% Plot the averages with no shifting
figure(5); clf(5); grid on; hold on;
xlim([0,27]);
ylim([-40,120]);
h1 = plot(xvec, Favg_nomag_mea, 'Color', colorsMat2(1,:), 'LineWidth',1,'LineStyle','--');
fill([xvec fliplr(xvec)],[Favg_nomag_mea + std_nomag_mea, fliplr(Favg_nomag_mea-std_nomag_mea)],...
    'b','FaceAlpha',0.2,'LineStyle','none');

h2 = plot(xvec, Favg_mag, 'Color', colorsMat2(2,:), 'LineWidth',2,'LineStyle','-');
fill([xvec fliplr(xvec)],[Favg_mag + std_mag, fliplr(Favg_mag-std_mag)],...
    'r','FaceAlpha',0.2,'LineStyle','none');

% h3 = plot(xvec, Favg_nomag_ea, 'Color', colorsMat2(3,:), 'LineWidth',3,'LineStyle',':');
% fill([xvec fliplr(xvec)],[Favg_nomag_ea + std_nomag_ea, fliplr(Favg_nomag_ea-std_nomag_ea)],...
%     'g','FaceAlpha',0.2,'LineStyle','none');

h4 = plot(xvec,-(Favg_nomag_mea-Favg_mag),'Color',colorsMat2(4,:),'LineWidth',4,'LineStyle',':');

% Plot differences
% h5 = plot(xvec,(Favg_nomag_mea-Favg_nomag_ea),'Color',colorsMat2(5,:),'LineWidth',4,'LineStyle',':');
% h6 = plot(xvec,-(Favg_nomag_ea-Favg_mag),'Color',colorsMat2(6,:),'LineWidth',4,'LineStyle',':');
legend([h1,h2,h4],'No Magnetic Guidance N=4','Magnetic Guidance N=4','Average Difference N=4');
% legend([h1,h2,h3,h4,h5,h6],{'Favg_{nomag_mea}','Favg_{mag}',...
%     'Favg_{nomag_ea}','Favgdiff_{nomagMEA2mag}','Favgdiff_{nomagMEA2nomagEA}',...
%     'Favgdiff_{nomagEA2mag}'});
xlabel('Insertion Depth [mm]'); ylabel('Average ||Force|| [mN]');