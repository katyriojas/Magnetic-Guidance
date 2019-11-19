% Katy Riojas and Trevor Bruns
% This script plots the averages from the robotic cadaver trials

% Last Updated: 11/19/19

%% Calculate the maximum X and Y of the plot (this will be the minimum
% linear displacement and minimum Fmag
min_robotic_cadaver_X = 1000;
max_robotic_cadaver_Y = 0;

for ii = 1:size(data_robotic_cadaver,2)
    min_robotic_cadaver_X = min([min_robotic_cadaver_X,...
                                 data_robotic_cadaver(ii).nomag.depth_insertion(robotic_cadaver_nomag_end(ii)),...
                                 data_robotic_cadaver(ii).mag.depth_insertion(robotic_cadaver_mag_end(ii))]);
    max_robotic_cadaver_Y = max([max_robotic_cadaver_Y,...
                                 data_robotic_cadaver(ii).nomag.Fmag_smooth(robotic_cadaver_nomag_end(ii)),...
                                 data_robotic_cadaver(ii).mag.Fmag_smooth(robotic_cadaver_mag_end(ii))]);
end

%% Interpolate so that we have 
xvec = linspace(0.02,min_robotic_cadaver_X,1000);

for ii = 1:size(data_robotic_cadaver,2)
    
    Fsmooth_cadaver_nomag(ii,:) = ...
        interp1(data_robotic_cadaver(ii).nomag.depth_insertion(1:robotic_cadaver_nomag_end(ii)),...
        data_robotic_cadaver(ii).nomag.Fmag_smooth(1:robotic_cadaver_nomag_end(ii)),...
        xvec);
    Fsmooth_cadaver_mag(ii,:) = ...
        interp1(data_robotic_cadaver(ii).mag.depth_insertion(1:robotic_cadaver_mag_end(ii)),...
        data_robotic_cadaver(ii).mag.Fmag_smooth(1:robotic_cadaver_mag_end(ii)),...
        xvec);
   
end


% Compute Averages and Standard Deviations
Favg_nomag = mean(Fsmooth_cadaver_nomag,1);
std_nomag = std(Fsmooth_cadaver_nomag);

Favg_mag = mean(Fsmooth_cadaver_mag,1);
std_mag = std(Fsmooth_cadaver_mag);

%% Plot the averages +/- one standard deviation
set(gca,'fontname','Arial');

figure(1); clf(1); grid on; hold on;
xlabel('Insertion Depth [mm]'); ylabel('Average ||Force|| [mN]');
xlim([0,min_robotic_cadaver_X]);
ylim([0,max_robotic_cadaver_Y]);

h1 = plot(xvec, Favg_nomag, 'Color', 'r', 'LineWidth',1,'LineStyle','--');
fill([xvec fliplr(xvec)],[Favg_nomag + std_nomag, fliplr(Favg_nomag-std_nomag)],...
    'r','FaceAlpha',0.2,'LineStyle','none');
h2 = plot(xvec, Favg_mag, 'Color', 'b', 'LineWidth',1,'LineStyle','-');
fill([xvec fliplr(xvec)],[Favg_mag + std_mag, fliplr(Favg_mag-std_mag)],...
    'b','FaceAlpha',0.2,'LineStyle','none');
legend([h1,h2],{'Robotic','Robotic & Magnetic Steering'});