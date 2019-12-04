% Katy Riojas and Trevor Bruns
% This script plots the averages from the robotic cadaver trials

% Last Updated: 11/19/19
regenerate_robotic_cadaver_data = false;

if regenerate_robotic_cadaver_data
    LoadRALData_Robotic_Cadaver;
elseif ~exist('data_robotic_cadaver','var') % if not already loaded
    load('data\cadaver\data_robotic_cadaver.mat'); % load already generated
end

%% Calculate the maximum X and Y of the plot (this will be the minimum
% linear displacement and minimum Fmag
min_robotic_cadaver_X = 1000;
max_robotic_cadaver_Y = 0;

for ii = 1:size(data_robotic_cadaver,2)
    min_robotic_cadaver_X = min([min_robotic_cadaver_X,...
                                 data_robotic_cadaver(ii).nomag_depth_insertion_trimmed(end),...
                                 data_robotic_cadaver(ii).mag_depth_insertion_trimmed(end)]);
    max_robotic_cadaver_Y = max([max_robotic_cadaver_Y,...
                                 max(data_robotic_cadaver(ii).nomag_Fmagsmooth_trimmed),...
                                 max(data_robotic_cadaver(ii).mag_Fmagsmooth_trimmed)]);
end

%% Interpolate so that we have 
xvec = linspace(0.02,min_robotic_cadaver_X,1000);

for ii = 1:size(data_robotic_cadaver,2)
    
    Fmag_cadaver_ug(ii,:) = ...
        interp1(data_robotic_cadaver(ii).nomag_depth_insertion_trimmed,...
                data_robotic_cadaver(ii).nomag_Fmagsmooth_trimmed,...
                xvec);
    Fmag_cadaver_g(ii,:) = ...
        interp1(data_robotic_cadaver(ii).mag_depth_insertion_trimmed,...
                data_robotic_cadaver(ii).mag_Fmagsmooth_trimmed,...
                xvec);
   
end


% Compute Averages and Standard Deviations
Favg_nomag = mean(Fmag_cadaver_ug,1);
std_nomag = std(Fmag_cadaver_ug);

Favg_mag = mean(Fmag_cadaver_g,1);
std_mag = std(Fmag_cadaver_g);

%% Plot the averages +/- one standard deviation
set(gca,'fontname','Arial');

figure(1); clf(1); grid on; hold on;
xlabel('Insertion Depth (mm)','FontWeight','bold'); 
ylabel('Average Force (mN)','FontWeight','bold');
xlim([0,min_robotic_cadaver_X]);
ylim([0,max_robotic_cadaver_Y]);

h1 = plot(xvec, Favg_nomag, 'Color', 'b', 'LineWidth',1,'LineStyle','--');
fill([xvec fliplr(xvec)],[Favg_nomag + std_nomag, fliplr(Favg_nomag-std_nomag)],...
    'b','FaceAlpha',0.2,'LineStyle','none');
h2 = plot(xvec, Favg_mag, 'Color', 'g', 'LineWidth',1,'LineStyle','-');
fill([xvec fliplr(xvec)],[Favg_mag + std_mag, fliplr(Favg_mag-std_mag)],...
    'g','FaceAlpha',0.2,'LineStyle','none');
legend([h1,h2],{'Robotic','Robotic & Magnetic Steering'},'Location','northwest');

set(0,'DefaultAxesFontSize',10);
set(0,'DefaultAxesFontName','Arial');
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 3.6 2];
saveas(fig,'saved figures\AvgCadavervs.LID.pdf');