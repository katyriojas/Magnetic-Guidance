% Data Analysis for Robotic Cadaver Trials
% Outputs:
    % Fig. 1: Plot Fmag vs. LID
    % Fig. 2: Plot Fx,Fy,Fz vs. LID
    % User can specify whether to regenerate phantom data, and whether to
    % update saved data_robotic_phantom.mat.

% Trevor Bruns and Katy Riojas
% Last Updated: 11/21/19

%% Regenerate data if requested
regenerate_robotic_cadaver_data = false;

if regenerate_robotic_cadaver_data
    LoadRALData_Robotic_Cadaver;
elseif ~exist('data_robotic_cadaver','var') % if not already loaded
    load('data\cadaver\data_robotic_cadaver.mat'); % load already generated
end

%% Initializations
addpath('functions');
alpha = 1; % line opacity
colorsMat = distinguishable_colors(6); % pull in a set of distinguishable colors
max_robotic_cadaver_Y = 0;
xyzColors = [1,0,0;0,1,0;0,0,1;...
            1,0,1;0,0.2,0.13;0,1,1];
line_width_smooth = 1.5;
line_width_raw = 0.75;

% Calc max y all
for ii = 1:size(data_robotic_cadaver,2)
    max_robotic_cadaver_Y = max([max_robotic_cadaver_Y,...
                            max(data_robotic_cadaver(ii).nomag.Fmag),...
                            max(data_robotic_cadaver(ii).mag.Fmag)]);
end

%% Plot 1: Fmag vs. LID Smoothed on top of Raw and plot trim points as scatter points
figure(10); clf(10);
sgtitle('Raw and Smoothed Fmag Cadaver Trials');
% First plot the unguided data

for ii = 1:size(data_robotic_cadaver,2)
    subplot(2,size(data_robotic_cadaver,2),ii); 
    grid on; hold on; ylim([0 max_robotic_cadaver_Y]);
    title(strcat('Trial ',num2str(ii)));
    if ii == 1
        ylabel('Unguided ||Force|| [mN]');
    end
    plot(data_robotic_cadaver(ii).nomag.depth_insertion, ...
         data_robotic_cadaver(ii).nomag.Fmag,...
        'Color', [colorsMat(1,:),  0.3*alpha], 'LineWidth',line_width_raw);
    plot(data_robotic_cadaver(ii).nomag.depth_insertion, ...
        data_robotic_cadaver(ii).nomag.Fmag_smooth,...
        'Color', colorsMat(1,:), 'LineWidth',line_width_smooth);
    scatter(data_robotic_cadaver(ii).nomag_depth_insertion_trimmed(end),...
            data_robotic_cadaver(ii).nomag_Fmag_trimmed(end),...
            10,'filled','r');
        
    % Plot the guided data on the bottom row
    subplot(2,size(data_robotic_cadaver,2),ii+size(data_robotic_cadaver,2)); 
    grid on; hold on; ylim([0 max_robotic_cadaver_Y]);
    title(strcat('Trial ',num2str(ii)));
    
    if ii == 1
       ylabel('Guided ||Force|| [mN]');
    end
    
    plot(data_robotic_cadaver(ii).mag.depth_insertion, ...
         data_robotic_cadaver(ii).mag.Fmag,...
        'Color', [colorsMat(1,:), 0.3*alpha], 'LineWidth',line_width_raw);
    plot(data_robotic_cadaver(ii).mag.depth_insertion, ...
        data_robotic_cadaver(ii).mag.Fmag_smooth,...
        'Color', colorsMat(1,:), 'LineWidth',line_width_smooth);
    scatter(data_robotic_cadaver(ii).mag_depth_insertion_trimmed(end),...
            data_robotic_cadaver(ii).mag_Fmag_trimmed(end),...
            10,'filled','r');   
    
end

%% Fig. 2: Plot Fx,Fy,Fz vs. LID
figure(2); clf(2); 
hold on; grid on;

for ii = 1:size(data_robotic_cadaver,2)
   
    % Plot unguided XYZ data               
                     plot(data_robotic_cadaver(ii).nomag.depth_insertion, data_robotic_cadaver(ii).nomag.Fx,...
                          'Color', xyzColors(1,:), 'LineWidth', line_width_raw);
    h(ii).nomag(1) = plot(data_robotic_cadaver(ii).nomag.depth_insertion, data_robotic_cadaver(ii).nomag.Fx_smooth,...
                          'Color', xyzColors(1,:), 'LineWidth', line_width_smooth);
                     plot(data_robotic_cadaver(ii).nomag.depth_insertion, data_robotic_cadaver(ii).nomag.Fy,...
                          'Color', xyzColors(2,:), 'LineWidth', line_width_raw);
    h(ii).nomag(2) = plot(data_robotic_cadaver(ii).nomag.depth_insertion, data_robotic_cadaver(ii).nomag.Fy_smooth,...
                          'Color', xyzColors(2,:), 'LineWidth', line_width_smooth);
                     plot(data_robotic_cadaver(ii).nomag.depth_insertion, data_robotic_cadaver(ii).nomag.Fz,...
                          'Color', xyzColors(3,:), 'LineWidth', line_width_raw);
    h(ii).nomag(3) = plot(data_robotic_cadaver(ii).nomag.depth_insertion, data_robotic_cadaver(ii).nomag.Fz_smooth,...
                          'Color', xyzColors(3,:), 'LineWidth', line_width_smooth);
    
    % Plot XYZ Guided Data
                   plot(data_robotic_cadaver(ii).mag.depth_insertion, data_robotic_cadaver(ii).mag.Fx,...
                        'Color', xyzColors(4,:), 'LineWidth', line_width_raw);
    h(ii).mag(1) = plot(data_robotic_cadaver(ii).mag.depth_insertion, data_robotic_cadaver(ii).mag.Fx_smooth,...
                        'Color', xyzColors(4,:), 'LineWidth', line_width_smooth);
                   plot(data_robotic_cadaver(ii).mag.depth_insertion, data_robotic_cadaver(ii).mag.Fy,...
                        'Color', xyzColors(5,:), 'LineWidth', line_width_raw);
    h(ii).mag(2) = plot(data_robotic_cadaver(ii).mag.depth_insertion, data_robotic_cadaver(ii).mag.Fy_smooth,...
                        'Color', xyzColors(5,:), 'LineWidth', line_width_smooth);
                   plot(data_robotic_cadaver(ii).mag.depth_insertion, data_robotic_cadaver(ii).mag.Fz,...
                        'Color', xyzColors(6,:), 'LineWidth', line_width_raw);
    h(ii).mag(3) = plot(data_robotic_cadaver(ii).mag.depth_insertion, data_robotic_cadaver(ii).mag.Fz_smooth,...
                        'Color', xyzColors(6,:), 'LineWidth', line_width_smooth);
end

title('Fx, Fy, and Fz vs Linear Insertion Depth')
xlabel('Linear Insertion Depth (mm)')
ylabel('Force (mN)')
legend([h(1).nomag,h(1).mag], [{'UG-Fx','UG-FY','UG-Fz'},{'G-Fx', 'G-Fy', 'G-Fz'}], 'Location', 'sw')