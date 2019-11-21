%% Compare binned phantom forces

% Output is a binned plot and we will save the binned data for use with the
% colorbar plot

% Trevor Bruns and Katy Riojas
% Last Updated: 11/19/19

% TODO: Talk about whether to use smoothed or raw data (the plot generated
% here looks slightly different because this uses the raw data)

regenerate_manual_data = false;
regenerate_phantom_data = false;
update_saved_phantom_structs = true; % default to true - need binned data for colorbar plot

if regenerate_manual_data
    LoadRALData_Manual; % regen data
elseif ~exist('data_manual_phantom','var')
    load('data\phantom\data_manual_phantom.mat'); % load already generated
end

if regenerate_phantom_data
    LoadRALData_Robotic_Phantom % regen data
elseif ~exist('data_robotic_phantom','var') % if not already loaded
    load('data\phantom\data_robotic_phantom.mat'); % load already generated
end
%% Initialize Variables for Binning and Plotting Data
% First compute the largest angular insertion depth to bound our bins
max_binnedX = 0;
for ii = 1:size(data_manual_phantom,2)
    max_binnedX = max([max_binnedX,...
                       max(data_manual_phantom(ii).interp_angdepth),...
                       max(data_robotic_phantom(ii).nomag_mea_interp_angdepth),...
                       max(data_robotic_phantom(ii).mag_interp_angdepth)]);
end

degspan = 1;
stepsize = 1; 
numsteps = ceil(max_binnedX/degspan);           
AOIvec = linspace(0,max_binnedX,numsteps); % this is our vector of interest

%% Intialize our vectors to be used later
% Outputs from this loop should be both the binned individual trials, as
% well as the average binned data
manual_binnedFmag = zeros(length(AOIvec),4);
robotic_ug_binnedFmag = zeros(length(AOIvec),4);
robotic_g_binnedFmag = zeros(length(AOIvec),4);

for ii = 2:length(AOIvec)
    for jj = 1:size(data_manual_phantom,2)
        % Find the manual data points
        idxman = find((data_manual_phantom(jj).interp_angdepth<AOIvec(ii))&...
                      (data_manual_phantom(jj).interp_angdepth>AOIvec(ii-1)));
        Fmag_manual_jj = data_manual_phantom(jj).Fmag_trimmed(idxman);
        manual_binnedFmag(ii,jj) = mean(Fmag_manual_jj);
        
        % Find phantom UG data points in bins
        idxug = find((data_robotic_phantom(jj).nomag_mea_interp_angdepth<AOIvec(ii))&...
                     (data_robotic_phantom(jj).nomag_mea_interp_angdepth>AOIvec(ii-1)));
        Fmag_robotic_nomag_jj = data_robotic_phantom(jj).nomag_mea.Fmag(idxug);
        robotic_ug_binnedFmag(ii,jj) = mean(Fmag_robotic_nomag_jj);
        
        % Find phantom G data points in bins
        idxg = find((data_robotic_phantom(jj).mag_interp_angdepth<AOIvec(ii))&...
                     (data_robotic_phantom(jj).mag_interp_angdepth>AOIvec(ii-1)));
        Fmag_robotic_mag_jj = data_robotic_phantom(jj).mag.Fmag(idxg);
        robotic_g_binnedFmag(ii,jj) = mean(Fmag_robotic_mag_jj);
    end
end

% Save data and check for nans
for jj = 1:size(data_manual_phantom,2)
    binned_man_jj = [AOIvec',manual_binnedFmag(:,jj)];
    binned_man_jj(any(isnan(binned_man_jj),2),:) = [];
    data_manual_phantom(jj).binned.angle = binned_man_jj(:,1);
    data_manual_phantom(jj).binned.Fmag = binned_man_jj(:,2);
    
    binned_ug_jj = [AOIvec',robotic_ug_binnedFmag(:,jj)];
    binned_ug_jj(any(isnan(binned_ug_jj),2),:) = [];
    data_robotic_phantom(jj).nomag_mea_binned.angle = binned_ug_jj(:,1);
    data_robotic_phantom(jj).nomag_mea_binned.Fmag = binned_ug_jj(:,2);
    
    binned_g_jj = [AOIvec',robotic_g_binnedFmag(:,jj)];
    binned_g_jj(any(isnan(binned_g_jj),2),:) = [];
    data_robotic_phantom(jj).mag_binned.angle = binned_g_jj(:,1);
    data_robotic_phantom(jj).mag_binned.Fmag = binned_g_jj(:,2);
end

% Find Binned means
manual_binnedFmag_avg     =  nanmean(manual_binnedFmag,2);
robotic_ug_binnedFmag_avg =  nanmean(robotic_ug_binnedFmag,2);
robotic_g_binnedFmag_avg  =  nanmean(robotic_g_binnedFmag,2);

%% Check Data Fits
figure(2); clf(2);
sgtitle('Binned Phantom Data Check');
subplot(1,3,1); grid on; hold on;
title('Manual Insertions');
for ii = 1:size(data_manual_phantom,2)
    plot(data_manual_phantom(ii).interp_angdepth,...
         data_manual_phantom(ii).Fmag_trimmed,...
         'Color','b','LineWidth',1,'LineStyle',':');
    plot(AOIvec-degspan/2, manual_binnedFmag_avg,...
         'Color', 'k','LineWidth',1,'LineStyle','-');
end

subplot(1,3,2); grid on; hold on;
title('UG Robotic Insertions');
for ii = 1:size(data_robotic_phantom,2)
    plot(data_robotic_phantom(ii).nomag_mea_interp_angdepth,...
         data_robotic_phantom(ii).nomag_mea.Fmag,...
         'Color','b','LineWidth',1,'LineStyle',':');
    plot(AOIvec-degspan/2, robotic_ug_binnedFmag_avg,...
         'Color', 'k','LineWidth',1,'LineStyle','-');
end

subplot(1,3,3); grid on; hold on;
title('G Robotic Insertions');
for ii = 1:size(data_robotic_phantom,2)
    plot(data_robotic_phantom(ii).mag_interp_angdepth,...
         data_robotic_phantom(ii).mag.Fmag,...
         'Color','b','LineWidth',1,'LineStyle',':');
    plot(AOIvec-degspan/2, robotic_g_binnedFmag_avg,...
         'Color', 'k','LineWidth',1,'LineStyle','-');
end

%% Plot average data comparison
figure(1); clf(1); grid on; hold on; xlim([0 400]);
xlabel('Angular Insertion Depth (\circ)'); ylabel('Average ||Force|| (mN)');
plot(AOIvec-degspan/2, manual_binnedFmag_avg, 'Color', 'r','LineWidth',1,'LineStyle','-');
plot(AOIvec-degspan/2, robotic_ug_binnedFmag_avg, 'Color', 'b','LineWidth',1,'LineStyle','-');
plot(AOIvec-degspan/2, robotic_g_binnedFmag_avg, 'Color', 'g','LineWidth',1,'LineStyle','-');
legend('Manual','Robotic','Robotic & Magnetic Steering');

%% Update Saved Phantom Data if it is called for
if update_saved_phantom_structs
   save('data\phantom\data_manual_phantom.mat','data_manual_phantom');
   save('data\phantom\data_robotic_phantom.mat','data_robotic_phantom');
end