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


%% Combine binned force measurements for each set of trials and compute statistics

clear phantom_stats;

phantom_stats.Fmag.bins = data_robotic_phantom(1).mag_binned.bins;
for i_bin = 1:length(phantom_stats.Fmag.bins)
    
    % create vectors containing all force measurements within the corresponding bin
    nomag_mea_binned(i_bin).Fmags = [];
    mag_binned(i_bin).Fmags = [];
    for i_trial = 1:length(data_robotic_phantom)
        % append forces from current trial
        nomag_mea_binned(i_bin).Fmags = [nomag_mea_binned(i_bin).Fmags; ...
                                      data_robotic_phantom(i_trial).nomag_mea.Fmag(data_robotic_phantom(i_trial).nomag_mea_binned.ind == i_bin)];

        mag_binned(i_bin).Fmags = [mag_binned(i_bin).Fmags; ...
                                data_robotic_phantom(i_trial).mag.Fmag(data_robotic_phantom(i_trial).mag_binned.ind == i_bin)];
    end

    manual_binned(i_bin).Fmags = [];
    for i_trial = 1:length(data_manual_phantom)
        % append forces from current trial
        manual_binned(i_bin).Fmags = [manual_binned(i_bin).Fmags; ...
                                      data_manual_phantom(i_trial).Fmag_trimmed(data_manual_phantom(i_trial).binned.ind == i_bin)];
    end
    
    % compute mean for current bin
    phantom_stats.Fmag.mean.nomag(i_bin)  = mean(nomag_mea_binned(i_bin).Fmags);
    phantom_stats.Fmag.mean.mag(i_bin)    = mean(mag_binned(i_bin).Fmags);
    phantom_stats.Fmag.mean.manual(i_bin) = mean(manual_binned(i_bin).Fmags);

    % compute standard deviation
    phantom_stats.Fmag.std.nomag(i_bin)  = std(nomag_mea_binned(i_bin).Fmags);
    phantom_stats.Fmag.std.mag(i_bin)    = std(mag_binned(i_bin).Fmags);
    phantom_stats.Fmag.std.manual(i_bin) = std(manual_binned(i_bin).Fmags);

    % perform t-test to determine if mean force with magnet is less than without
%     if (~isnan(phantom_stats.Fmag.mean.mag(i_bin)) && ~isnan(phantom_stats.Fmag.mean.nomag(i_bin)) ) % must have values for both
        [phantom_stats.Fmag.diff.h(i_bin),  phantom_stats.Fmag.diff.p(i_bin), phantom_stats.Fmag.diff.ci(i_bin,:)] = ...
            ttest2( mag_binned(i_bin).Fmags, nomag_mea_binned(i_bin).Fmags, 'Tail','left', 'Vartype','unequal');

        % compute mean difference between manual/robotic
        phantom_stats.Fmag.diff.mean(i_bin) = phantom_stats.Fmag.mean.mag(i_bin) - phantom_stats.Fmag.mean.nomag(i_bin);
%     end

end

% remove NaNs and convert to logicals
phantom_stats.Fmag.diff.h(isnan(phantom_stats.Fmag.diff.h)) = 0;
phantom_stats.Fmag.diff.h = logical(phantom_stats.Fmag.diff.h);

%% Check Data Fits
% figure(2); clf(2);
% sgtitle('Binned Phantom Data Check');
% subplot(1,3,1); grid on; hold on;
% title('Manual Insertions');
% for i_bin = 1:size(data_manual_phantom,2)
%     plot(data_manual_phantom(i_bin).interp_angdepth,...
%          data_manual_phantom(i_bin).Fmag_trimmed,...
%          'Color','b','LineWidth',1,'LineStyle',':');
%     plot(AOIvec-degspan/2, manual_binnedFmag_avg,...
%          'Color', 'k','LineWidth',1,'LineStyle','-');
% end
% 
% subplot(1,3,2); grid on; hold on;
% title('UG Robotic Insertions');
% for i_bin = 1:size(data_robotic_phantom,2)
%     plot(data_robotic_phantom(i_bin).nomag_mea_interp_angdepth,...
%          data_robotic_phantom(i_bin).nomag_mea.Fmag,...
%          'Color','b','LineWidth',1,'LineStyle',':');
%     plot(AOIvec-degspan/2, robotic_ug_binnedFmag_avg,...
%          'Color', 'k','LineWidth',1,'LineStyle','-');
% end
% 
% subplot(1,3,3); grid on; hold on;
% title('G Robotic Insertions');
% for i_bin = 1:size(data_robotic_phantom,2)
%     plot(data_robotic_phantom(i_bin).mag_interp_angdepth,...
%          data_robotic_phantom(i_bin).mag.Fmag,...
%          'Color','b','LineWidth',1,'LineStyle',':');
%     plot(AOIvec-degspan/2, robotic_g_binnedFmag_avg,...
%          'Color', 'k','LineWidth',1,'LineStyle','-');
% end

%% Plot average data comparison
if exist('hf_avg_binned','var')
    if isvalid(hf_avg_binned)
        close(hf_avg_binned)
    end
end

hf_avg_binned = figure;
hf_avg_binned.WindowState = "maximized";

h_ax(1) = subplot_er(2,1,1);
grid on; hold on; %xlim([0 400]);

% plot means
plot(phantom_stats.Fmag.bins, phantom_stats.Fmag.mean.manual, 'Color', 'r','LineWidth',1.5);
plot(phantom_stats.Fmag.bins, phantom_stats.Fmag.mean.nomag,  'Color', 'b','LineWidth',1.5);
plot(phantom_stats.Fmag.bins, phantom_stats.Fmag.mean.mag,    'Color', 'g','LineWidth',1.5);

% plot standard deviations as shaded regions around means
range = find(~isnan(phantom_stats.Fmag.mean.manual),1) : (length(phantom_stats.Fmag.bins) - find(~isnan( fliplr(phantom_stats.Fmag.mean.manual)),1)); 
fill([phantom_stats.Fmag.bins(range), fliplr(phantom_stats.Fmag.bins(range))],...
     [phantom_stats.Fmag.mean.manual(range) + phantom_stats.Fmag.std.manual(range), fliplr(phantom_stats.Fmag.mean.manual(range) - phantom_stats.Fmag.std.manual(range))],...
     'r', 'FaceAlpha',0.15, 'EdgeColor','none');

range = find(~isnan(phantom_stats.Fmag.mean.nomag),1) : (length(phantom_stats.Fmag.bins) - find(~isnan( fliplr(phantom_stats.Fmag.mean.nomag)),1)); 
fill([phantom_stats.Fmag.bins(range), fliplr(phantom_stats.Fmag.bins(range))],...
     [phantom_stats.Fmag.mean.nomag(range) + phantom_stats.Fmag.std.nomag(range), fliplr(phantom_stats.Fmag.mean.nomag(range) - phantom_stats.Fmag.std.nomag(range))],...
     'b', 'FaceAlpha',0.15, 'EdgeColor','none');

range = find(~isnan(phantom_stats.Fmag.mean.mag),1) : (length(phantom_stats.Fmag.bins) - find(~isnan( fliplr(phantom_stats.Fmag.mean.mag)),1)); 
fill([phantom_stats.Fmag.bins(range), fliplr(phantom_stats.Fmag.bins(range))],...
     [phantom_stats.Fmag.mean.mag(range) + phantom_stats.Fmag.std.mag(range), fliplr(phantom_stats.Fmag.mean.mag(range) - phantom_stats.Fmag.std.mag(range))],...
     'g', 'FaceAlpha',0.15, 'EdgeColor','none');

% mark trial end depths
for i_trial = 1:length(data_robotic_phantom)
    last_ind = data_manual_phantom(i_trial).binned.ind(end);
    scatter(phantom_stats.Fmag.bins(last_ind), phantom_stats.Fmag.mean.manual(last_ind), 100, 'r', 'd', 'filled');

    last_ind = data_robotic_phantom(i_trial).nomag_mea_binned.ind(end);
    scatter(phantom_stats.Fmag.bins(last_ind), phantom_stats.Fmag.mean.nomag(last_ind),  100, 'b', 'd', 'filled');

    last_ind = data_robotic_phantom(i_trial).mag_binned.ind(end-1); % TODO: fix NaN in interp_angdepth
    scatter(phantom_stats.Fmag.bins(last_ind), phantom_stats.Fmag.mean.mag(last_ind),    100, 'g', 'd', 'filled');
end

ylabel('||F|| (mN)');
legend('Manual','Robotic','Robotic & Magnetic Steering', 'Location','nw');



h_ax(2) = subplot_er(2,1,2);
grid on; hold on;
xlabel('Angular Insertion Depth (\circ)'); 
ylabel('\Delta ||F|| (mN)');
not_nan = ~isnan(phantom_stats.Fmag.diff.mean);
plot(phantom_stats.Fmag.bins(not_nan), phantom_stats.Fmag.diff.mean(not_nan), 'Color', 'k','LineWidth',1.5);
area(phantom_stats.Fmag.bins(not_nan), phantom_stats.Fmag.diff.mean(not_nan), 'FaceColor',[0,0,0], 'FaceAlpha',0.05, 'EdgeColor','none');
scatter(phantom_stats.Fmag.bins(phantom_stats.Fmag.diff.h), phantom_stats.Fmag.diff.mean(phantom_stats.Fmag.diff.h), 100, 'm', '*')


linkaxes(h_ax, 'x');
xlim([0, phantom_stats.Fmag.bins(end)])

%% Update Saved Phantom Data if it is called for
% if update_saved_phantom_structs
%    save('data\phantom\data_manual_phantom.mat','data_manual_phantom');
%    save('data\phantom\data_robotic_phantom.mat','data_robotic_phantom');
% end