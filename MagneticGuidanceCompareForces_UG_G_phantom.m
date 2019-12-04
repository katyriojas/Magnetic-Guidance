% Data Analysis for Robotic Phantom Insertions
% Outputs:
    % Fig. 1: Plot Fmag vs. LID
    % Fig. 2: Plot Fx,Fy,Fz vs. LID
    % Fig. 3: Plot AID vs. LID
    % Fig. 4: Plot Fmag vs. LID
    % Fig. 5: Plot Fx,Fy,Fz vs. AID
    % User can specify whether to regenerate phantom data, and whether to
    % update saved data_robotic_phantom.mat.

% Trevor Bruns and Katy Riojas
% Last Updated: 11/21/19

%% Regenerate Phantom data if requested
regenerate_phantom_data = false;

% guided/unguided phantom data
if regenerate_phantom_data
    LoadRALData_Robotic_Phantom % regen data
elseif ~exist('data_robotic_phantom','var') % if not already loaded
    load('data\phantom\data_robotic_phantom.mat'); % load already generated
end

%% Initializations
addpath('functions');
colors = distinguishable_colors(2*length(data_robotic_phantom)+1);
alpha = 1; % reduce transparency of unguided plot lines
xyzColors = [1,0,0;0,1,0;0,0,1;...
            1,0,1;0,0.2,0.13;0,1,1];
line_width_smooth = 1.5;
line_width_raw = 0.75;

%% Plot -> Fmag vs Linear Insertion Depth
% Figure: force magnitude
figure(1); clf(1); hold on; grid on;

for ii = 1:size(data_robotic_phantom,2)

                    plot(data_robotic_phantom(ii).nomag_ea.depth_insertion, data_robotic_phantom(ii).nomag_ea.Fmag,...
                        'Color',[colors(ii,:), 0.3*alpha],'LineWidth',line_width_raw);
  h(ii).nomag_ea =  plot(data_robotic_phantom(ii).nomag_ea.depth_insertion, data_robotic_phantom(ii).nomag_ea.Fmag_smooth,...
                        'Color', colors(ii,:),'LineStyle',':', 'LineWidth',line_width_smooth);
 
                    plot(data_robotic_phantom(ii).nomag_mea.depth_insertion, data_robotic_phantom(ii).nomag_mea.Fmag,...
                        'Color', [colors(ii,:),0.3*alpha],'LineWidth',line_width_raw);
  h(ii).nomag_mea = plot(data_robotic_phantom(ii).nomag_mea.depth_insertion, data_robotic_phantom(ii).nomag_mea.Fmag_smooth,...
                        'Color',colors(ii,:),'LineStyle','--','LineWidth',line_width_smooth);

                    plot(data_robotic_phantom(ii).mag.depth_insertion, data_robotic_phantom(ii).mag.Fmag,...
                        'Color',[colors(ii,:),0.3*alpha],'LineWidth',line_width_raw);
  h(ii).mag =       plot(data_robotic_phantom(ii).mag.depth_insertion, data_robotic_phantom(ii).mag.Fmag_smooth,...
                        'Color', colors(ii,:),'LineWidth',line_width_smooth);

end

title('Guided vs Unguided Insertion')
xlabel('Insertion Depth [mm]')
ylabel('||Force|| [mN]')

legend([h.nomag_ea, h.nomag_mea, h.mag],{'nomag1-ea', 'nomag2-ea','nomag3-ea','nomag4-ea',...
                                         'nomag1-mea','nomag2-mea','nomag3-mea','nomag4-mea',...
                                         'mag1','mag2','mag3','mag4'});
                                     
%% Plot -> XYZ Force vs. Linear Insertion Depth
figure(2); clf(2); 
hold on; grid on;

for ii = 1:size(data_robotic_phantom,2)
   
    % Plot unguided XYZ data               
                     plot(data_robotic_phantom(ii).nomag_mea.depth_insertion, data_robotic_phantom(ii).nomag_mea.Fx,...
                          'Color', xyzColors(1,:), 'LineWidth', line_width_raw);
    h(ii).nomag(1) = plot(data_robotic_phantom(ii).nomag_mea.depth_insertion, data_robotic_phantom(ii).nomag_mea.Fx_smooth,...
                          'Color', xyzColors(1,:), 'LineWidth', line_width_smooth);
                     plot(data_robotic_phantom(ii).nomag_mea.depth_insertion, data_robotic_phantom(ii).nomag_mea.Fy,...
                          'Color', xyzColors(2,:), 'LineWidth', line_width_raw);
    h(ii).nomag(2) = plot(data_robotic_phantom(ii).nomag_mea.depth_insertion, data_robotic_phantom(ii).nomag_mea.Fy_smooth,...
                          'Color', xyzColors(2,:), 'LineWidth', line_width_smooth);
                     plot(data_robotic_phantom(ii).nomag_mea.depth_insertion, data_robotic_phantom(ii).nomag_mea.Fz,...
                          'Color', xyzColors(3,:), 'LineWidth', line_width_raw);
    h(ii).nomag(3) = plot(data_robotic_phantom(ii).nomag_mea.depth_insertion, data_robotic_phantom(ii).nomag_mea.Fz_smooth,...
                          'Color', xyzColors(3,:), 'LineWidth', line_width_smooth);
    
    % Plot XYZ Guided Data
                   plot(data_robotic_phantom(ii).mag.depth_insertion, data_robotic_phantom(ii).mag.Fx,...
                        'Color', xyzColors(4,:), 'LineWidth', line_width_raw);
    h(ii).mag(1) = plot(data_robotic_phantom(ii).mag.depth_insertion, data_robotic_phantom(ii).mag.Fx_smooth,...
                        'Color', xyzColors(4,:), 'LineWidth', line_width_smooth);
                   plot(data_robotic_phantom(ii).mag.depth_insertion, data_robotic_phantom(ii).mag.Fy,...
                        'Color', xyzColors(5,:), 'LineWidth', line_width_raw);
    h(ii).mag(2) = plot(data_robotic_phantom(ii).mag.depth_insertion, data_robotic_phantom(ii).mag.Fy_smooth,...
                        'Color', xyzColors(5,:), 'LineWidth', line_width_smooth);
                   plot(data_robotic_phantom(ii).mag.depth_insertion, data_robotic_phantom(ii).mag.Fz,...
                        'Color', xyzColors(6,:), 'LineWidth', line_width_raw);
    h(ii).mag(3) = plot(data_robotic_phantom(ii).mag.depth_insertion, data_robotic_phantom(ii).mag.Fz_smooth,...
                        'Color', xyzColors(6,:), 'LineWidth', line_width_smooth);
end

title('Fx, Fy, and Fz vs Linear Insertion Depth')
xlabel('Linear Insertion Depth (mm)')
ylabel('Force (mN)')
legend([h(1).nomag,h(1).mag,], [{'UG-Fx','UG-FY','UG-Fz'},{'G-Fx', 'G-Fy', 'G-Fz'}], 'Location', 'sw')

%% Plot the AID vs. LID to make sure that seemed reasonable.
figure(3); clf(3); hold on; grid on;

for ii = 1:size(data_robotic_phantom,2)
    
    h_nomag(ii) = plot(data_robotic_phantom(ii).nomag_mea.depth_insertion,...
                       data_robotic_phantom(ii).nomag_mea_interp_angdepth,...
                       'Color',colors(ii*2,:), 'LineStyle',':', 'LineWidth', line_width_smooth);
    h_mag(ii)   = plot(data_robotic_phantom(ii).mag.depth_insertion,...
                       data_robotic_phantom(ii).mag_interp_angdepth,...
                       'Color',colors(ii*2-1,:), 'LineWidth', line_width_smooth);

end

xlabel('Actuator Insertion Distance (mm)')
ylabel('Angular Insertion Depth (deg)')

clear labels;
labels.mag = {'mag1','mag2','mag3','mag4'};
labels.nomag = {'nomag1','nomag2','nomag3','nomag4'};
legend([h_nomag(1:ii),h_mag(1:ii)], [labels.nomag(1:ii),labels.mag(1:ii)], 'Location','nw')

%% Plot Fmag vs. AID
figure(4); clf(4);
hold on; grid on;

for ii=1:size(data_robotic_phantom,2)
   
    % Plot Unguided Fmag
                plot(data_robotic_phantom(ii).nomag_mea_interp_angdepth, data_robotic_phantom(ii).nomag_mea.Fmag,...
                'Color', [colors(ii,:), 0.3], 'LineWidth', 1, 'LineStyle',':', 'LineWidth', line_width_raw);
    h_nomag(ii) = plot(data_robotic_phantom(ii).nomag_mea_interp_angdepth, data_robotic_phantom(ii).nomag_mea.Fmag_smooth,...
                'Color', [colors(ii,:), 1],   'LineWidth', 1, 'LineStyle',':', 'LineWidth', line_width_smooth);
    
    % Plot guided Fmag
                plot(data_robotic_phantom(ii).mag_interp_angdepth, data_robotic_phantom(ii).mag.Fmag,...
                'Color', [colors(ii,:), 0.3], 'LineWidth', line_width_raw);
    h_mag(ii) = plot(data_robotic_phantom(ii).mag_interp_angdepth, data_robotic_phantom(ii).mag.Fmag_smooth,...
                'Color', [colors(ii,:), 1], 'LineWidth', line_width_smooth);
    
end

title('Force vs Angular Insertion Depth')
xlabel('Angular insertion depth (deg)')
ylabel('||Force|| (mN)')

clear labels;
labels.mag = {'mag1','mag2','mag3','mag4'};
labels.nomag = {'nomag1','nomag2','nomag3','nomag4'};
legend([h_nomag(1:ii),h_mag(1:ii)], [labels.nomag(1:ii),labels.mag(1:ii)], 'Location','nw')

%% Plot -> XYZ Force vs. Angular insertion depth
figure(5); clf(5); 
hold on; grid on;

for ii = 1:size(data_robotic_phantom,2)
   
    % Plot unguided XYZ data               
                     plot(data_robotic_phantom(ii).nomag_mea_interp_angdepth, data_robotic_phantom(ii).nomag_mea.Fx,...
                          'Color', xyzColors(1,:), 'LineWidth', line_width_raw);
    h(ii).nomag(1) = plot(data_robotic_phantom(ii).nomag_mea_interp_angdepth, data_robotic_phantom(ii).nomag_mea.Fx_smooth,...
                          'Color', xyzColors(1,:), 'LineWidth', line_width_smooth);
                     plot(data_robotic_phantom(ii).nomag_mea_interp_angdepth, data_robotic_phantom(ii).nomag_mea.Fy,...
                          'Color', xyzColors(2,:), 'LineWidth', line_width_raw);
    h(ii).nomag(2) = plot(data_robotic_phantom(ii).nomag_mea_interp_angdepth, data_robotic_phantom(ii).nomag_mea.Fy_smooth,...
                          'Color', xyzColors(2,:), 'LineWidth', line_width_smooth);
                     plot(data_robotic_phantom(ii).nomag_mea_interp_angdepth, data_robotic_phantom(ii).nomag_mea.Fz,...
                          'Color', xyzColors(3,:), 'LineWidth', line_width_raw);
    h(ii).nomag(3) = plot(data_robotic_phantom(ii).nomag_mea_interp_angdepth, data_robotic_phantom(ii).nomag_mea.Fz_smooth,...
                          'Color', xyzColors(3,:), 'LineWidth', line_width_smooth);
    
    % Plot XYZ Guided Data
                   plot(data_robotic_phantom(ii).mag_interp_angdepth, data_robotic_phantom(ii).mag.Fx,...
                        'Color', xyzColors(4,:), 'LineWidth', line_width_raw);
    h(ii).mag(1) = plot(data_robotic_phantom(ii).mag_interp_angdepth, data_robotic_phantom(ii).mag.Fx_smooth,...
                        'Color', xyzColors(4,:), 'LineWidth', line_width_smooth);
                   plot(data_robotic_phantom(ii).mag_interp_angdepth, data_robotic_phantom(ii).mag.Fy,...
                        'Color', xyzColors(5,:), 'LineWidth', line_width_raw);
    h(ii).mag(2) = plot(data_robotic_phantom(ii).mag_interp_angdepth, data_robotic_phantom(ii).mag.Fy_smooth,...
                        'Color', xyzColors(5,:), 'LineWidth', line_width_smooth);
                   plot(data_robotic_phantom(ii).mag_interp_angdepth, data_robotic_phantom(ii).mag.Fz,...
                        'Color', xyzColors(6,:), 'LineWidth', line_width_raw);
    h(ii).mag(3) = plot(data_robotic_phantom(ii).mag_interp_angdepth, data_robotic_phantom(ii).mag.Fz_smooth,...
                        'Color', xyzColors(6,:), 'LineWidth', line_width_smooth);
end

title('Force vs Linear Insertion Depth')
xlabel('Linear Insertion Depth (mm)')
ylabel('||Force|| (mN)')
legend([h(1).nomag,h(1).mag,], [{'UG-Fx','UG-FY','UG-Fz'},{'G-Fx', 'G-Fy', 'G-Fz'}], 'Location', 'sw')