%% Plot Results from Manual Insertion Trials
% This script generates the following plots:
    % 1) Raw Fmag Data vs. Time all Trials
    %   - note raw Fmag data is used for trimming
    % 2) Trimmed Fmag vs. Time Phantom
    % 3) Trimmed Fmag vs Time Cadaver
    % 4) xyz Fmag vs. Time Phantom and Cadaver

% Trevor Bruns and Katy Riojas
% Last Updated: 11/19/19


%% Regenerate data if requested
regenerate_manual_data = false;

if regenerate_manual_data
    LoadRALData_Manual; % regen data
elseif ~exist('data_manual_phantom','var')
    load('data\phantom\data_manual_phantom.mat'); % load already generated
end

%% Plotting Variables
addpath('functions');
alpha = 1; % line opacity
colorsRed = [1,0,0;...
             1,0,0.3;...
             0.76,0.23,0.13;...
              1,0,0]; % pull in a set of distinguishable colors
max_robotic_cadaver_Y = 0;
xyzColors = [1,0,0;0,1,0;0,0,1;...
            1,0,1;0,0.2,0.13;0,1,1];
line_width_smooth = 1.5;
line_width_raw = 0.75;
line_width_vline = 2;
line_style_mat = {'-','--',':'};

%% Calculate Max Y for both phantom and cadaver
max_manual_phantomY = 0;
max_manual_cadaverY = 0;
max_manual_phantomX = 0;
max_manual_cadaverX = 0;

% Use trimmed results for calculation
% Calc max y phantom
for ii = 1:length(filepaths_manual_phantom)
%     max_manual_phantomX = max(max_manual_phantomX,max(data_manual_phantom(ii).nano.time));
%     max_manual_phantomY = max(max_manual_phantomY,max(data_manual_phantom(ii).nano.Fmag));
    max_manual_phantomX = max(max_manual_phantomX,max(data_manual_phantom(ii).time_trimmed));
    max_manual_phantomY = max(max_manual_phantomY,max(data_manual_phantom(ii).Fmag_trimmed));
end

% Calc max y cadaver
for ii = 1:length(filepaths_manual_cadaver)
%     max_manual_cadaverX = max(max_manual_cadaverX,max(data_manual_cadaver(ii).nano.time));
%     max_manual_cadaverY = max(max_manual_cadaverY,max(data_manual_cadaver(ii).nano.Fmag));
    max_manual_cadaverX = max(max_manual_cadaverX,max(data_manual_cadaver(ii).time_trimmed));
    max_manual_cadaverY = max(max_manual_cadaverY,max(data_manual_cadaver(ii).Fmag_trimmed));
end

%% Plot 1: Force magnitude vs. Time for cropping
figure(1); clf(1);
sgtitle('Manual Data: Phantom Fmag vs. Time');

for ii = 1:size(data_manual_phantom,2)
    subplot(2,4,ii); grid on; hold on; xlabel('Time (s)'); 
    ylim([0,max_manual_phantomY]);
    title(strcat('Trial ',num2str(ii)));
    if ii == 1
        ylabel ('Fmag Phantom ||F|| (mN)');
    end
    
    % Plot raw (non-smoothed) Fmag vs. Time
    h1(ii) = plot(data_manual_phantom(ii).nano.time,...
                  data_manual_phantom(ii).nano.Fmag,'Color','b');
    
    % Plot trimmed start and end times
    h2(ii) = line([data_manual_phantom(ii).nano.time(data_manual_phantom(ii).trim_idx(1)),...
                   data_manual_phantom(ii).nano.time(data_manual_phantom(ii).trim_idx(1))],...
                  [0,max_manual_phantomY],...
                  'Color','k','LineStyle','--','LineWidth',line_width_vline);
             line([data_manual_phantom(ii).nano.time(data_manual_phantom(ii).trim_idx(end)),...
                   data_manual_phantom(ii).nano.time(data_manual_phantom(ii).trim_idx(end))],...
                   [0,max_manual_phantomY],...
                   'Color','k','LineStyle','--','LineWidth',line_width_vline);
    
    % Put scatter point at trim point location
    scatter(data_manual_phantom(ii).nano.time(data_manual_phantom(ii).trim_idx(1)),...
            data_manual_phantom(ii).Fmag_trimmed(1),...
            10,'filled','r');
    scatter(data_manual_phantom(ii).nano.time(data_manual_phantom(ii).trim_idx(end)),...
            data_manual_phantom(ii).Fmag_trimmed(end),...
            10,'filled','r');
    
    % Plot release time for each trial
    h3(ii) = line([releaseTimes(ii),releaseTimes(ii)],...
                  [0,max_manual_phantomY],...
                  'Color','g','LineStyle','--','LineWidth',line_width_vline);
    
    legend([h1(ii),h2(ii),h3(ii)],{'Fmag','Trim Times','Release Time'});
    
    if ii ~= size(data_manual_phantom,2) % then plot the cadaver data on bottom row
        subplot(2,3,ii+3); grid on; hold on; xlabel('Time (s)'); 
        ylim([0,max_manual_cadaverY]);
        title(strcat('Trial ',num2str(ii)));
        if ii == 1
            ylabel ('Fmag Cadaver ||F|| (mN)');
        end

        % Plot raw (non-smoothed) Fmag vs. Time
        h1(ii) = plot(data_manual_cadaver(ii).nano.time,...
                      data_manual_cadaver(ii).nano.Fmag,'Color','b');

        % Plot trimmed start and end times
        h2(ii) = line([data_manual_cadaver(ii).nano.time(data_manual_cadaver(ii).trim_idx(1)),...
                       data_manual_cadaver(ii).nano.time(data_manual_cadaver(ii).trim_idx(1))],...
                       [0,max_manual_cadaverY],...
                       'Color','k','LineStyle','--','LineWidth',line_width_vline);
                 line([data_manual_cadaver(ii).nano.time(data_manual_cadaver(ii).trim_idx(end)),...
                       data_manual_cadaver(ii).nano.time(data_manual_cadaver(ii).trim_idx(end))],...
                      [0,max_manual_cadaverY],...
                      'Color','k','LineStyle','--','LineWidth',line_width_vline);

        % Put scatter point at trim point location
        scatter(data_manual_cadaver(ii).nano.time(data_manual_cadaver(ii).trim_idx(1)),...
                data_manual_cadaver(ii).Fmag_trimmed(1),...
                10,'filled','r');
        scatter(data_manual_cadaver(ii).nano.time(data_manual_cadaver(ii).trim_idx(end)),...
                data_manual_cadaver(ii).Fmag_trimmed(end),...
                10,'filled','r');

        % Plot release time for each trial (only have relase time for first
        % cadaver for some reason
        %h3(ii) = line([releaseTimes(ii),releaseTimes(ii)],[0,max_manual_cadaverY],'Color','g','LineStyle','--','LineWidth',line_width_vline);

        legend([h1(ii),h2(ii),h3(ii)],{'Fmag','Trim Times','Release Time'});
    end
end

%% Plot 2: Plot the Force Magnitude vs. Time for Manual Phantom Trials
figure(2); clf(2); grid on; hold on; 
title('Trimmed Phantom Manual Trials N = 4');
xlabel('Insertion Time (s)'); ylabel('Force (mN)');
xlim([0,max_manual_phantomX]);
ylim([0,max_manual_phantomY]);

for ii = 1:size(data_manual_phantom,2)
    plot(data_manual_phantom(ii).time_trimmed,...
         data_manual_phantom(ii).Fmag_trimmed,...
         'Color',colorsRed(ii,:),'LineWidth',line_width_raw);
end

legend('Trial 1','Trial 2','Trial 3','Trial 4');

%% Plot the Force Magnitude vs. Time for Manual Cadaver Trials
figure(3); clf(3); grid on; hold on; 
% title('Trimmed Cadaver Manual Trials N = 3');
xlabel('Insertion Time (s)','FontWeight','bold'); 
ylabel('Force (mN)','FontWeight','bold');
xlim([0,max_manual_cadaverX]);
ylim([0,max_manual_cadaverY]);

for ii = 1:size(data_manual_cadaver,2)
    plot(data_manual_cadaver(ii).time_trimmed,...
         data_manual_cadaver(ii).Fmag_trimmed,...
         'Color',colorsRed(ii,:),'LineWidth',line_width_raw,'LineStyle',...
         line_style_mat{ii});
end

legend('Trial 1','Trial 2','Trial 3','Location','northwest');
set(0,'DefaultAxesFontSize',10);
set(0,'DefaultAxesFontName','Arial');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 3.6 2];
saveas(fig,'saved figures\ManualCadavervs.Time.pdf');

%% Plot 3: XYZ Plots using trimmed data
figure(4); clf(4); hold on; grid on;
sgtitle('Fx, Fy, Fz vs Time');

for ii = 1:size(data_manual_phantom,2)
   
    subplot(2,4,ii); hold on; grid on;
    title(strcat('Trial ',num2str(ii))); 
    xlabel('Time (s)'); ylabel('Phantom Force (mN)');  
    xlim([0,max_manual_phantomX]);
    ylim([0,60]); % hard coded - check this number if data changes
    
             plot(data_manual_phantom(ii).nano.time, data_manual_phantom(ii).nano.Fx,...
                  'Color', xyzColors(1,:), 'LineWidth', line_width_raw);
    hx(ii) = plot(data_manual_phantom(ii).nano.time, data_manual_phantom(ii).nano.Fx_smooth,...
                  'Color', xyzColors(1,:), 'LineWidth', line_width_smooth);
                     plot(data_manual_phantom(ii).nano.time, data_manual_phantom(ii).nano.Fy,...
                  'Color', xyzColors(2,:), 'LineWidth', line_width_raw);
    hy(ii) = plot(data_manual_phantom(ii).nano.time, data_manual_phantom(ii).nano.Fy_smooth,...
                  'Color', xyzColors(2,:), 'LineWidth', line_width_smooth);
                     plot(data_manual_phantom(ii).nano.time, data_manual_phantom(ii).nano.Fz,...
                  'Color', xyzColors(3,:), 'LineWidth', line_width_raw);
    hz(ii) = plot(data_manual_phantom(ii).nano.time, data_manual_phantom(ii).nano.Fz_smooth,...
                  'Color', xyzColors(3,:), 'LineWidth', line_width_smooth);
    legend([hx(ii),hy(ii),hz(ii)], {'Fx','Fy','Fz'}, 'Location', 'sw');
    
    if ii ~= size(data_manual_phantom,2)
        subplot(2,3,ii+3); hold on; grid on;
        title(strcat('Trial ',num2str(ii))); 
        xlabel('Time (s)'); ylabel('Cadaver Force (mN)');
        xlim([0,max_manual_cadaverX]);
        ylim([0,max_manual_cadaverY]);
                 plot(data_manual_cadaver(ii).nano.time, data_manual_cadaver(ii).nano.Fx,...
                      'Color', xyzColors(1,:), 'LineWidth', line_width_raw);
        hx(ii) = plot(data_manual_cadaver(ii).nano.time, data_manual_cadaver(ii).nano.Fx_smooth,...
                      'Color', xyzColors(1,:), 'LineWidth', line_width_smooth);
                 plot(data_manual_cadaver(ii).nano.time, data_manual_cadaver(ii).nano.Fy,...
                      'Color', xyzColors(2,:), 'LineWidth', line_width_raw);
        hy(ii) = plot(data_manual_cadaver(ii).nano.time, data_manual_cadaver(ii).nano.Fy_smooth,...
                      'Color', xyzColors(2,:), 'LineWidth', line_width_smooth);
                 plot(data_manual_cadaver(ii).nano.time, data_manual_cadaver(ii).nano.Fz,...
                      'Color', xyzColors(3,:), 'LineWidth', line_width_raw);
        hz(ii) = plot(data_manual_cadaver(ii).nano.time, data_manual_cadaver(ii).nano.Fz_smooth,...
                      'Color', xyzColors(3,:), 'LineWidth', line_width_smooth);
        legend([hx(ii),hy(ii),hz(ii)], {'Fx','Fy','Fz'}, 'Location', 'sw');
    end
end