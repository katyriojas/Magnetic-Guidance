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
    [phantom_stats.Fmag.diff.h(i_bin),  phantom_stats.Fmag.diff.p(i_bin), phantom_stats.Fmag.diff.ci(i_bin,:)] = ...
            ttest2( mag_binned(i_bin).Fmags, nomag_mea_binned(i_bin).Fmags, 'Tail','left', 'Vartype','unequal', 'Alpha', 0.01);

    % compute mean difference between manual/robotic
    phantom_stats.Fmag.diff.mean(i_bin) = phantom_stats.Fmag.mean.mag(i_bin) - phantom_stats.Fmag.mean.nomag(i_bin);

    % compute RMS standard deviation of the difference
    phantom_stats.Fmag.diff.std(i_bin) = sqrt( phantom_stats.Fmag.std.nomag(i_bin)^2 + phantom_stats.Fmag.std.mag(i_bin)^2);

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

line_width = 2;
alpha_std = 0.1;

h_ax(1) = subplot_er(2,1,1);
grid on; hold on;
set(gca,'FontSize',9,'FontName','Times');
xlim([0, phantom_stats.Fmag.bins(end)+10]);
% title(strcat(sprintf('Mean Forces in Phantom (Bin Size = %i', phantom_stats.Fmag.bins(2)-phantom_stats.Fmag.bins(1)), '\circ)'))

% plot means
plot(phantom_stats.Fmag.bins, phantom_stats.Fmag.mean.manual, 'Color', 'r','LineWidth',line_width);
plot(phantom_stats.Fmag.bins, phantom_stats.Fmag.mean.nomag,  'Color', 'b','LineWidth',line_width);
plot(phantom_stats.Fmag.bins, phantom_stats.Fmag.mean.mag,    'Color', 'g','LineWidth',line_width);

% plot standard deviations as shaded regions around means
range = find(~isnan(phantom_stats.Fmag.mean.manual),1) : (length(phantom_stats.Fmag.bins) - find(~isnan( fliplr(phantom_stats.Fmag.mean.manual)),1)); 
fill([phantom_stats.Fmag.bins(range), fliplr(phantom_stats.Fmag.bins(range))],...
     [phantom_stats.Fmag.mean.manual(range) + phantom_stats.Fmag.std.manual(range), fliplr(phantom_stats.Fmag.mean.manual(range) - phantom_stats.Fmag.std.manual(range))],...
     'r', 'FaceAlpha',alpha_std, 'EdgeColor','none');

range = find(~isnan(phantom_stats.Fmag.mean.nomag),1) : (length(phantom_stats.Fmag.bins) - find(~isnan( fliplr(phantom_stats.Fmag.mean.nomag)),1)); 
fill([phantom_stats.Fmag.bins(range), fliplr(phantom_stats.Fmag.bins(range))],...
     [phantom_stats.Fmag.mean.nomag(range) + phantom_stats.Fmag.std.nomag(range), fliplr(phantom_stats.Fmag.mean.nomag(range) - phantom_stats.Fmag.std.nomag(range))],...
     'b', 'FaceAlpha',alpha_std, 'EdgeColor','none');

range = find(~isnan(phantom_stats.Fmag.mean.mag),1) : (length(phantom_stats.Fmag.bins) - find(~isnan( fliplr(phantom_stats.Fmag.mean.mag)),1)); 
fill([phantom_stats.Fmag.bins(range), fliplr(phantom_stats.Fmag.bins(range))],...
     [phantom_stats.Fmag.mean.mag(range) + phantom_stats.Fmag.std.mag(range), fliplr(phantom_stats.Fmag.mean.mag(range) - phantom_stats.Fmag.std.mag(range))],...
     'g', 'FaceAlpha',alpha_std, 'EdgeColor','none');

 ms = 20;
% mark trial end depths
for i_trial = 1:length(data_robotic_phantom)
    last_ind = data_manual_phantom(i_trial).binned.ind(end);
    scatter(phantom_stats.Fmag.bins(last_ind), phantom_stats.Fmag.mean.manual(last_ind), ms, 'r', 'd', 'filled','MarkerEdgeColor','k');

    last_ind = data_robotic_phantom(i_trial).nomag_mea_binned.ind(end);
    scatter(phantom_stats.Fmag.bins(last_ind), phantom_stats.Fmag.mean.nomag(last_ind),  ms, 'b', 'd', 'filled','MarkerEdgeColor','k');

    last_ind = data_robotic_phantom(i_trial).mag_binned.ind(end-1); % TODO: fix NaN in interp_angdepth
    scatter(phantom_stats.Fmag.bins(last_ind), phantom_stats.Fmag.mean.mag(last_ind),    ms, 'g', 'd', 'filled','MarkerEdgeColor','k');
end

ylabel('||F|| (mN)','FontName','Times','FontSize',9,'FontWeight','bold');
legend('Manual','Robotic','Robotic & Magnetic Steering', 'Location','nw',...
    'FontName','Times','FontSize',9);
% ylim([-10, h_ax(1).YLim(2)])
% ylim([0 120])

% Plot delta F between guided/unguided robotic insertion means
h_ax(2) = subplot_er(2,1,2);
grid on; hold on;
set(gca,'FontSize',9,'FontName','Times');

% delta F
not_nan = ~isnan(phantom_stats.Fmag.diff.mean);
plot(phantom_stats.Fmag.bins(not_nan), phantom_stats.Fmag.diff.mean(not_nan), 'Color', 'k','LineWidth',line_width);

% standard deviation
% range = find(~isnan(phantom_stats.Fmag.mean.manual),1) : (length(phantom_stats.Fmag.bins) - find(~isnan( fliplr(phantom_stats.Fmag.mean.manual)),1)); 
fill([phantom_stats.Fmag.bins(not_nan), fliplr(phantom_stats.Fmag.bins(not_nan))],...
     [phantom_stats.Fmag.diff.mean(not_nan) + phantom_stats.Fmag.diff.std(not_nan), fliplr(phantom_stats.Fmag.diff.mean(not_nan) - phantom_stats.Fmag.diff.std(not_nan))],...
     'k', 'FaceAlpha',alpha_std, 'EdgeColor','none');

% area(phantom_stats.Fmag.bins(not_nan), phantom_stats.Fmag.diff.mean(not_nan), 'FaceColor',[0,0,0], 'FaceAlpha',0.05, 'EdgeColor','none');

% mark t-test significant points
scatter(phantom_stats.Fmag.bins(phantom_stats.Fmag.diff.h), phantom_stats.Fmag.diff.mean(phantom_stats.Fmag.diff.h), ms, 'm', '*')

xlabel('Angular Insertion Depth (\circ)','FontName','Times','FontSize',9,'FontWeight','bold'); 
ylabel('\Delta||F|| (mN)','FontName','Times','FontSize',9,'FontWeight','bold');

linkaxes(h_ax, 'x');
xlim([0, phantom_stats.Fmag.bins(end)+10])
ylim([min(phantom_stats.Fmag.diff.mean(not_nan) - phantom_stats.Fmag.diff.std(not_nan) - 10), max(phantom_stats.Fmag.diff.mean(not_nan) + phantom_stats.Fmag.diff.std(not_nan) + 10)]) 
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperSize = [8.5 11];
fig.PaperPosition = [0 0 3.5 2.5];
saveas(fig,'saved figures\PhantomBinned_Comparison.pdf');

%% Plot p-values from t-test
% figure(20); clf(20);
% bar(phantom_stats.Fmag.bins, phantom_stats.Fmag.diff.p)


%% Update Saved Phantom Data if it is called for
% if update_saved_phantom_structs
%    save('data\phantom\data_manual_phantom.mat','data_manual_phantom');
%    save('data\phantom\data_robotic_phantom.mat','data_robotic_phantom');
% end