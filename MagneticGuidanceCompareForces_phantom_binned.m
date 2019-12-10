%% Compare binned phantom forces

% Output is a binned plot and we will save the binned data for use with the
% colorbar plot

% Trevor Bruns and Katy Riojas
% Last Updated: December 2019

regenerate_manual_data = false;
regenerate_phantom_data = false;

%% Initializations

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

% RALData_Binning;

colors = distinguishable_colors(2*length(data_robotic_phantom)+1);
alpha = 1; % reduce transparency of unguided plot lines
xyzColors = [1,0,0;0,1,0;0,0,1;...
            1,0,1;0,0.2,0.13;0,1,1];
line_width_smooth = 1;
line_width_raw = 0.75;


%% Plot Fmag vs. AID
figure(4); clf(4);
hold on; grid on;

for ii=1:length(data_robotic_phantom)
    
%%% Unguided  %%%

    % Plot Unguided Fmag
                  plot(data_robotic_phantom(ii).nomag_mea_interp_angdepth, data_robotic_phantom(ii).nomag_mea.Fmag,...
                  'Color', [colors(ii,:), 0.3], 'LineStyle',':', 'LineWidth', line_width_raw);
    h_nomag(ii) = plot(data_robotic_phantom(ii).nomag_mea_interp_angdepth, data_robotic_phantom(ii).nomag_mea.Fmag_smooth,...
                  'Color', [colors(ii,:), 1],   'LineStyle',':', 'LineWidth', line_width_smooth);

    % find NaN indices (i.e. angles without force measurements)
    nan_inds = isnan(data_robotic_phantom(ii).nomag_mea_binned.Fmean);

    % plot Fmean for each bin   
    plot(data_robotic_phantom(ii).nomag_mea_binned.bins(~nan_inds), data_robotic_phantom(ii).nomag_mea_binned.Fmean(~nan_inds),...
        'Color', [colors(ii,:), 1],   'LineStyle',':', 'LineWidth', 1.5);
    
    % mark cutoff points
    if ~isnan(trim.nomag(ii).ind_cutoff)
        scatter(data_robotic_phantom(ii).nomag_mea_interp_angdepth(trim.nomag(ii).ind_cutoff),...
                data_robotic_phantom(ii).nomag_mea.Fmag_smooth(trim.nomag(ii).ind_cutoff), 50, colors(ii,:), 'd', 'filled');
    end

%%% Guided %%%

    % Plot guided Fmag
                plot(data_robotic_phantom(ii).mag_interp_angdepth, data_robotic_phantom(ii).mag.Fmag,...
                'Color', [colors(ii,:), 0.3], 'LineWidth', line_width_raw);
    h_mag(ii) = plot(data_robotic_phantom(ii).mag_interp_angdepth, data_robotic_phantom(ii).mag.Fmag_smooth,...
                'Color', [colors(ii,:), 1], 'LineWidth', line_width_smooth);

    % find NaN indices (i.e. angles without force measurements)
    nan_inds = isnan(data_robotic_phantom(ii).mag_binned.Fmean);

    % plot Fmean for each bin
    plot(data_robotic_phantom(ii).mag_binned.bins(~nan_inds), data_robotic_phantom(ii).mag_binned.Fmean(~nan_inds),...
        'Color', [colors(ii,:), 1], 'LineWidth', 1.5); 

    % mark cutoff (trim) points
    if ~isnan(trim.mag(ii).ind_cutoff)
        scatter(data_robotic_phantom(ii).mag_interp_angdepth(trim.mag(ii).ind_cutoff),...
                data_robotic_phantom(ii).mag.Fmag_smooth(trim.mag(ii).ind_cutoff), 50, colors(ii,:), 'd', 'filled');
    end
end

title('Force vs Angular Insertion Depth')
xlabel('Angular insertion depth (deg)')
ylabel('||Force|| (mN)')

clear labels;
labels.mag = {'mag1','mag2','mag3','mag4'};
labels.nomag = {'nomag1','nomag2','nomag3','nomag4'};
legend([h_nomag(1:ii),h_mag(1:ii)], [labels.nomag(1:ii),labels.mag(1:ii)], 'Location','nw')



%% Plot Averaged Fmag vs AID (trimmed)

if exist('hf_avg_binned_trimmed','var')
    if isvalid(hf_avg_binned_trimmed)
        close(hf_avg_binned_trimmed)
    end
end

hf_avg_binned_trimmed = figure;
% hf_avg_binned_trimmed.WindowState = "maximized";

line_width = 2;
alpha_std = 0.22;
ms = 100; % default 20;
ms2 = 30; % default 10;
fs = 28; %pnt

h_ax_t(1) = subplot_er(2,1,1);
grid on; hold on;
set(gca,'FontSize',fs,'FontName','Times');
title(strcat(sprintf('Mean Forces in Phantom - Trimmed (Bin Size = %i', phantom_stats_trim.Fmag.bins(2)-phantom_stats_trim.Fmag.bins(1)), '\circ)'))
set(gca,'xticklabel',{[]});

% plot means
plot(phantom_stats_trim.Fmag.bins, phantom_stats_trim.Fmag.mean.manual, 'Color', 'r','LineWidth',line_width);
plot(phantom_stats_trim.Fmag.bins, phantom_stats_trim.Fmag.mean.nomag,  'Color', 'b','LineWidth',line_width);
plot(phantom_stats_trim.Fmag.bins, phantom_stats_trim.Fmag.mean.mag,    'Color', 'g','LineWidth',line_width);

% plot standard deviations as shaded regions around means
range = find(~isnan(phantom_stats_trim.Fmag.mean.manual),1) : (length(phantom_stats_trim.Fmag.bins) - find(~isnan( fliplr(phantom_stats_trim.Fmag.mean.manual)),1)); 
fill([phantom_stats_trim.Fmag.bins(range), fliplr(phantom_stats_trim.Fmag.bins(range))],...
     [phantom_stats_trim.Fmag.mean.manual(range) + phantom_stats_trim.Fmag.std.manual(range), fliplr(phantom_stats_trim.Fmag.mean.manual(range) - phantom_stats_trim.Fmag.std.manual(range))],...
     'r', 'FaceAlpha',alpha_std, 'EdgeColor','none');

range = find(~isnan(phantom_stats_trim.Fmag.mean.nomag),1) : (length(phantom_stats_trim.Fmag.bins) - find(~isnan( fliplr(phantom_stats_trim.Fmag.mean.nomag)),1)); 
fill([phantom_stats_trim.Fmag.bins(range), fliplr(phantom_stats_trim.Fmag.bins(range))],...
     [phantom_stats_trim.Fmag.mean.nomag(range) + phantom_stats_trim.Fmag.std.nomag(range), fliplr(phantom_stats_trim.Fmag.mean.nomag(range) - phantom_stats_trim.Fmag.std.nomag(range))],...
     'b', 'FaceAlpha',alpha_std, 'EdgeColor','none');

range = find(~isnan(phantom_stats_trim.Fmag.mean.mag),1) : (length(phantom_stats_trim.Fmag.bins) - find(~isnan( fliplr(phantom_stats_trim.Fmag.mean.mag)),1)); 
fill([phantom_stats_trim.Fmag.bins(range), fliplr(phantom_stats_trim.Fmag.bins(range))],...
     [phantom_stats_trim.Fmag.mean.mag(range) + phantom_stats_trim.Fmag.std.mag(range), fliplr(phantom_stats_trim.Fmag.mean.mag(range) - phantom_stats_trim.Fmag.std.mag(range))],...
     'g', 'FaceAlpha',alpha_std, 'EdgeColor','none');

% mark trial end depths
for i_trial = 1:length(data_robotic_phantom_trim)
    last_ind = data_manual_phantom(i_trial).binned.ind(end);
    scatter(phantom_stats_trim.Fmag.bins(last_ind), phantom_stats_trim.Fmag.mean.manual(last_ind), ms, 'r', 'd', 'filled','MarkerEdgeColor','k');

    last_ind = data_robotic_phantom_trim(i_trial).nomag_mea_binned.ind(end);
    scatter(phantom_stats_trim.Fmag.bins(last_ind), phantom_stats_trim.Fmag.mean.nomag(last_ind),  ms, 'b', 'd', 'filled','MarkerEdgeColor','k');

    last_ind = data_robotic_phantom_trim(i_trial).mag_binned.ind(end);
    scatter(phantom_stats_trim.Fmag.bins(last_ind), phantom_stats_trim.Fmag.mean.mag(last_ind),    120, 'g', 'd', 'filled');
end

ylabel('||F|| (mN)','FontName','Times','FontWeight','bold','FontSize',fs);
legend('Manual','Robotic','Robotic &\newlineMagnetic Steering', 'Location','nw','FontSize',fs,'FontName','Times');
ylim([-10, h_ax_t(1).YLim(2)])
yticks([0,25,50,75,100]);

% Plot delta F between guided/unguided robotic insertion means
h_ax_t(2) = subplot_er(2,1,2);
grid on; hold on;
set(gca,'FontSize',fs,'FontName','Times');

% delta F
not_nan = ~isnan(phantom_stats_trim.Fmag.diff.mean);
plot(phantom_stats_trim.Fmag.bins(not_nan), phantom_stats_trim.Fmag.diff.mean(not_nan), 'Color', 'k','LineWidth',line_width);

% standard deviation
fill([phantom_stats_trim.Fmag.bins(not_nan), fliplr(phantom_stats_trim.Fmag.bins(not_nan))],...
     [phantom_stats_trim.Fmag.diff.mean(not_nan) + phantom_stats_trim.Fmag.diff.std(not_nan), fliplr(phantom_stats_trim.Fmag.diff.mean(not_nan) - phantom_stats_trim.Fmag.diff.std(not_nan))],...
     'k', 'FaceAlpha',alpha_std, 'EdgeColor','none');

% plot t-test confidence interval as shaded regions around means
% fill([phantom_stats_trim.Fmag.bins(not_nan), fliplr(phantom_stats_trim.Fmag.bins(not_nan))],...
%      [phantom_stats_trim.Fmag.diff.ci(not_nan,1)', fliplr(phantom_stats_trim.Fmag.diff.ci(not_nan,2)')],...
%      'k', 'FaceAlpha',alpha_std, 'EdgeColor','none');

% mark t-test significant points
H = scatter(phantom_stats_trim.Fmag.bins(phantom_stats_trim.Fmag.diff.h), phantom_stats_trim.Fmag.diff.mean(phantom_stats_trim.Fmag.diff.h), ms2, 'm', 'o', 'LineWidth',0.5);
legend(H,'Reject Null Hypothesis', 'Location','nw','FontSize',fs,'FontName','Times');

xlabel('Angular Insertion Depth (\circ)','FontWeight','bold','FontSize',fs,'FontName','Times'); 
ylabel('\Delta||F|| (mN)','FontWeight','bold','FontSize',fs,'FontName','Times');

linkaxes(h_ax_t, 'x');
xlim([0, phantom_stats_trim.Fmag.bins(end)+5])
% xlim([0,400])
ylim([min(phantom_stats_trim.Fmag.diff.mean(not_nan) - phantom_stats_trim.Fmag.diff.std(not_nan) - 2), max(phantom_stats_trim.Fmag.diff.mean(not_nan) + phantom_stats_trim.Fmag.diff.std(not_nan) + 2)]) 
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperSize = [8.5 11];
fig.PaperPosition = [0 0 10.5 10.25]; % 10.5, 10.25 for poster; 3.5, 2.5 for paper
saveas(fig,'saved figures\binned_phantom.pdf');


%% Plot Averaged Fmag vs AID (untrimmed)

if exist('hf_avg_binned','var')
    if isvalid(hf_avg_binned)
        close(hf_avg_binned)
    end
end

hf_avg_binned = figure;
% hf_avg_binned.WindowState = "maximized";

line_width = 2;
alpha_std = 0.1;

h_ax(1) = subplot_er(2,1,1);
grid on; hold on;
set(gca,'FontSize',9,'FontName','Times');
title(strcat(sprintf('Mean Forces in Phantom - Untrimmed (Bin Size = %i', phantom_stats.Fmag.bins(2)-phantom_stats.Fmag.bins(1)), '\circ)'))
xlim([0, phantom_stats.Fmag.bins(end)+10]);

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

%% Fmag vs Normalized Insertion Depth (all trials)

% figure(18); clf(18);
% hold on; grid on;
% title('Phantom: ||F|| vs Normalized Insertion Depth');
% 
% line_width = 1;
% 
% for i_trial = 1:length(data_manual_phantom)
%     % manual
%     plot([data_manual_phantom(i_trial).normbin_Fmag.pre.mean; data_manual_phantom(i_trial).normbin_Fmag.post.mean],...
%         'Color', 'r', 'LineWidth', line_width);
%     % unguided
%     plot([data_robotic_phantom(i_trial).nomag_normbin_Fmag.pre.mean; data_robotic_phantom(i_trial).nomag_normbin_Fmag.post.mean],...
%         'Color', 'b', 'LineWidth', line_width);
%     % guided
%     plot([data_robotic_phantom(i_trial).mag_normbin_Fmag.pre.mean; data_robotic_phantom(i_trial).mag_normbin_Fmag.post.mean],...
%         'Color', 'g', 'LineWidth', line_width);
% end
% 
% % mark basal turn point
% xline(n_pre_bins, '--k', 'Basal Turn', 'LineWidth',2, 'LabelHorizontalAlignment','center', 'LabelVerticalAlignment','middle');
% 
% xlabel('Normalized Angular Insertion Depth');
% ylabel('||F|| (mN)');
% set(gca,'xticklabel',{[]})
% 
% legend('Manual','Robotic','Robotic & Magnetic Steering', 'Location','nw');



%% Fmag vs Normalized Insertion Depth (averages)

% if exist('hf_normbin','var')
%     if isvalid(hf_normbin)
%         close(hf_normbin)
%     end
% end
% 
% hf_normbin = figure;
% hf_normbin.WindowState = "maximized";
% 
% line_width = 2;
% alpha_std = 0.15;
% 
% % h_ax(1) = subplot_er(2,1,1);
% h_ax(1) = subplot(2,1,1);
% grid on; hold on; %xlim([0 400]);
% title(sprintf('Phantom: ||F|| vs Normalized Insertion Depth (pre/post basal turn bins = [%i, %i])', n_pre_bins, n_post_bins));
% 
% 
% % plot means
% plot(phantom_stats.Fmag.normbin.mean.manual, 'Color', 'r','LineWidth',line_width);
% plot(phantom_stats.Fmag.normbin.mean.nomag,  'Color', 'b','LineWidth',line_width);
% plot(phantom_stats.Fmag.normbin.mean.mag,    'Color', 'g','LineWidth',line_width);
% 
% % plot standard deviations as shaded regions around means
% n_normbins = n_pre_bins + n_post_bins;
% X = 1:n_normbins;
% 
% not_nan = ~isnan(phantom_stats.Fmag.normbin.mean.manual);
% fill([X(not_nan), fliplr(X(not_nan))],...
%      [phantom_stats.Fmag.normbin.mean.manual(not_nan) + phantom_stats.Fmag.normbin.std.manual(not_nan), fliplr(phantom_stats.Fmag.normbin.mean.manual(not_nan) - phantom_stats.Fmag.normbin.std.manual(not_nan))],...
%      'r', 'FaceAlpha',alpha_std, 'EdgeColor','none');
% 
% not_nan = ~isnan(phantom_stats.Fmag.normbin.mean.nomag);
% fill([X(not_nan), fliplr(X(not_nan))],...
%      [phantom_stats.Fmag.normbin.mean.nomag(not_nan) + phantom_stats.Fmag.normbin.std.nomag(not_nan), fliplr(phantom_stats.Fmag.normbin.mean.nomag(not_nan) - phantom_stats.Fmag.normbin.std.nomag(not_nan))],...
%      'b', 'FaceAlpha',alpha_std, 'EdgeColor','none');
% 
% not_nan = ~isnan(phantom_stats.Fmag.normbin.mean.mag);
% fill([X(not_nan), fliplr(X(not_nan))],...
%      [phantom_stats.Fmag.normbin.mean.mag(not_nan) + phantom_stats.Fmag.normbin.std.mag(not_nan), fliplr(phantom_stats.Fmag.normbin.mean.mag(not_nan) - phantom_stats.Fmag.normbin.std.mag(not_nan))],...
%      'g', 'FaceAlpha',alpha_std, 'EdgeColor','none');
% 
% 
% ylabel('||F|| (mN)');
% % legend('Robotic','Robotic & Magnetic Steering', 'Location','nw');
% legend('Manual','Robotic','Robotic & Magnetic Steering', 'Location','nw');
% ylim([-10, h_ax(1).YLim(2)])
% 
% % mark basal turn point
% xline(n_pre_bins, '--k', 'Basal Turn', 'LineWidth',2, 'LabelHorizontalAlignment','center', 'LabelVerticalAlignment','middle');
% set(gca,'xticklabel',{[]})
% 
% 
% 
% % Plot delta F between guided/unguided robotic insertion means
% h_ax(2) = subplot(2,1,2);
% % h_ax(2) = subplot_er(2,1,2);
% grid on; hold on;
% 
% % delta F
% plot(phantom_stats.Fmag.normbin.diff.mean, 'Color', 'k','LineWidth',line_width);
% 
% % standard deviation
% not_nan = ~isnan(phantom_stats.Fmag.normbin.diff.mean);
% fill([X(not_nan), fliplr(X(not_nan))],...
%      [phantom_stats.Fmag.normbin.diff.mean(not_nan) + phantom_stats.Fmag.normbin.diff.std(not_nan), fliplr(phantom_stats.Fmag.normbin.diff.mean(not_nan) - phantom_stats.Fmag.normbin.diff.std(not_nan))],...
%      'k', 'FaceAlpha',alpha_std, 'EdgeColor','none');
% 
% % mark t-test significant points
% scatter(X(phantom_stats.Fmag.normbin.diff.h), phantom_stats.Fmag.normbin.diff.mean(phantom_stats.Fmag.normbin.diff.h), 40, 'm', 'o', 'LineWidth',1.7)
% 
% % mark basal turn point
% xline(n_pre_bins, '--k', 'Basal Turn', 'LineWidth',2, 'LabelHorizontalAlignment','center', 'LabelVerticalAlignment','middle');
% 
% xlabel('Normalized Angular Insertion Depth');
% ylabel('\Delta ||F|| (mN)');
% set(gca,'xticklabel',{[]})
% 
% linkaxes(h_ax, 'x');


%% Plot p-values from t-test
% figure(20); clf(20);
% bar(phantom_stats_trim.Fmag.bins, phantom_stats_trim.Fmag.diff.p)
