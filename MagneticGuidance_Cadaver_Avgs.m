% Katy Riojas and Trevor Bruns
% This script plots the averages from the robotic cadaver trials

% Last Updated: 11/19/19
regenerate_robotic_cadaver_data = true;

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

%% Binning - Linear Depth (trimmed)

linspan = 0.125; % [mm] width of each bin

% compute the largest linear insertion depth to bound our bins
max_depth = 0;
for ii = 1:length(data_robotic_cadaver)
    max_depth = max([max_depth,...
                     max(data_robotic_cadaver(ii).nomag_depth_insertion_trimmed),...
                     max(data_robotic_cadaver(ii).mag_depth_insertion_trimmed)]);
end
num_bins = ceil(max_depth/linspan);

% create vector of the bin edges
bin_edges = 0:linspan:num_bins*linspan;

% create vector of the bin centers
bins = bin_edges(2:end) - linspan/2;

% sort into bins
for ii = 1:length(data_robotic_cadaver)

    % determine the bin for each interp_angdepth point
    [~,~,data_robotic_cadaver(ii).nomag_binned.ind] = histcounts(data_robotic_cadaver(ii).nomag_depth_insertion_trimmed, bin_edges);
    [~,~,data_robotic_cadaver(ii).mag_binned.ind]   = histcounts(data_robotic_cadaver(ii).mag_depth_insertion_trimmed, bin_edges);

    % compute mean of all the force measurements within each bin (NaN for bins with no measurements)
    data_robotic_cadaver(ii).nomag_binned.Fmean = zeros(size(bins));
    data_robotic_cadaver(ii).mag_binned.Fmean   = zeros(size(bins));
    for jj = 1:length(bins)
        data_robotic_cadaver(ii).nomag_binned.Fmean(jj) = ...
            mean( data_robotic_cadaver(ii).nomag_Fmag_trimmed( data_robotic_cadaver(ii).nomag_binned.ind == jj ) );

        data_robotic_cadaver(ii).mag_binned.Fmean(jj) = ...
            mean( data_robotic_cadaver(ii).mag_Fmag_trimmed( data_robotic_cadaver(ii).mag_binned.ind == jj ) );
    end

    % also save the actual bins and edges
    data_robotic_cadaver(ii).nomag_binned.bins  = bins;
    data_robotic_cadaver(ii).nomag_binned.edges = bin_edges;
    data_robotic_cadaver(ii).mag_binned.bins  = bins;
    data_robotic_cadaver(ii).mag_binned.edges = bin_edges;
end


%% Combine binned force measurements for each set of trials and compute statistics

clear cadaver_stats;

cadaver_stats.Fmag.bins = data_robotic_cadaver(1).mag_binned.bins;

for i_bin = 1:length(cadaver_stats.Fmag.bins)
    
    % create vectors containing all force measurements within the corresponding bin
    nomag_binned(i_bin).Fmags = [];
      mag_binned(i_bin).Fmags = [];

    for i_trial = 1:length(data_robotic_cadaver)
        % append forces from current trial
        nomag_binned(i_bin).Fmags = [nomag_binned(i_bin).Fmags; ...
                                     data_robotic_cadaver(i_trial).nomag_Fmag_trimmed( data_robotic_cadaver(i_trial).nomag_binned.ind == i_bin )];

        mag_binned(i_bin).Fmags = [mag_binned(i_bin).Fmags; ...
                                   data_robotic_cadaver(i_trial).mag_Fmag_trimmed( data_robotic_cadaver(i_trial).mag_binned.ind == i_bin )];
    end
    
    % compute mean for current bin
    cadaver_stats.Fmag.mean.nomag(i_bin)  = mean(nomag_binned(i_bin).Fmags);
    cadaver_stats.Fmag.mean.mag(i_bin)    = mean(mag_binned(i_bin).Fmags);

    % compute standard deviation
    cadaver_stats.Fmag.std.nomag(i_bin)  = std(nomag_binned(i_bin).Fmags);
    cadaver_stats.Fmag.std.mag(i_bin)    = std(mag_binned(i_bin).Fmags);

    % perform t-test to determine if mean force with magnet is less than without
    [cadaver_stats.Fmag.diff.h(i_bin),  cadaver_stats.Fmag.diff.p(i_bin), cadaver_stats.Fmag.diff.ci(i_bin,:)] = ...
            ttest2( mag_binned(i_bin).Fmags, nomag_binned(i_bin).Fmags, 'Tail','left', 'Vartype','unequal', 'Alpha', 0.05);

    % compute mean difference between manual/robotic
    cadaver_stats.Fmag.diff.mean(i_bin) = cadaver_stats.Fmag.mean.mag(i_bin) - cadaver_stats.Fmag.mean.nomag(i_bin);

    % compute standard error of the difference between means
    cadaver_stats.Fmag.diff.std(i_bin) = sqrt( (cadaver_stats.Fmag.std.nomag(i_bin)^2) / length(nomag_binned(i_bin).Fmags) ...
                                               + (cadaver_stats.Fmag.std.mag(i_bin)^2) / length(mag_binned(i_bin).Fmags) );
end

% remove NaNs and convert to logicals
cadaver_stats.Fmag.diff.h(isnan(cadaver_stats.Fmag.diff.h)) = 0;
cadaver_stats.Fmag.diff.h = logical(cadaver_stats.Fmag.diff.h);


%% Plot Averaged Fmag vs AID (trimmed)

if exist('hf_avg_cadaver_binned','var')
    if isvalid(hf_avg_cadaver_binned)
        close(hf_avg_cadaver_binned)
    end
end

hf_avg_cadaver_binned = figure;
% hf_avg_cadaver_binned.WindowState = "maximized";

line_width = 2;
alpha_std = 0.22;

h_ax_t(1) = subplot_er(2,1,1);
grid on; hold on;
title(strcat(sprintf('Mean Forces in Cadaver - Trimmed (Bin Size = %.3f', cadaver_stats.Fmag.bins(2)-cadaver_stats.Fmag.bins(1)), 'mm)'))

% plot means
plot(cadaver_stats.Fmag.bins, cadaver_stats.Fmag.mean.nomag,  'Color', 'b','LineWidth',line_width);
plot(cadaver_stats.Fmag.bins, cadaver_stats.Fmag.mean.mag,    'Color', 'g','LineWidth',line_width);

% plot standard deviations as shaded regions around means
range = find(~isnan(cadaver_stats.Fmag.mean.nomag),1) : (length(cadaver_stats.Fmag.bins) - find(~isnan( fliplr(cadaver_stats.Fmag.mean.nomag)),1)); 
fill([cadaver_stats.Fmag.bins(range), fliplr(cadaver_stats.Fmag.bins(range))],...
     [cadaver_stats.Fmag.mean.nomag(range) + cadaver_stats.Fmag.std.nomag(range), fliplr(cadaver_stats.Fmag.mean.nomag(range) - cadaver_stats.Fmag.std.nomag(range))],...
     'b', 'FaceAlpha',alpha_std, 'EdgeColor','none');

range = find(~isnan(cadaver_stats.Fmag.mean.mag),1) : (length(cadaver_stats.Fmag.bins) - find(~isnan( fliplr(cadaver_stats.Fmag.mean.mag)),1)); 
fill([cadaver_stats.Fmag.bins(range), fliplr(cadaver_stats.Fmag.bins(range))],...
     [cadaver_stats.Fmag.mean.mag(range) + cadaver_stats.Fmag.std.mag(range), fliplr(cadaver_stats.Fmag.mean.mag(range) - cadaver_stats.Fmag.std.mag(range))],...
     'g', 'FaceAlpha',alpha_std, 'EdgeColor','none');

% mark trial end depths
for i_trial = 1:length(data_robotic_cadaver)
    last_ind = data_robotic_cadaver(i_trial).nomag_binned.ind(end);
    scatter(cadaver_stats.Fmag.bins(last_ind), cadaver_stats.Fmag.mean.nomag(last_ind),  120, 'b', 'd', 'filled');

    last_ind = data_robotic_cadaver(i_trial).mag_binned.ind(end);
    scatter(cadaver_stats.Fmag.bins(last_ind), cadaver_stats.Fmag.mean.mag(last_ind),    120, 'g', 'd', 'filled');
end

ylabel('||F|| (mN)');
legend('Robotic','Robotic + Magnetic Steering', 'Location','nw');
ylim([-10, h_ax_t(1).YLim(2)])

% Plot delta F between guided/unguided robotic insertion means
h_ax_t(2) = subplot_er(2,1,2);
grid on; hold on;

% delta F
not_nan = ~isnan(cadaver_stats.Fmag.diff.mean);
plot(cadaver_stats.Fmag.bins(not_nan), cadaver_stats.Fmag.diff.mean(not_nan), 'Color', 'k','LineWidth',line_width);

% standard deviation
fill([cadaver_stats.Fmag.bins(not_nan), fliplr(cadaver_stats.Fmag.bins(not_nan))],...
     [cadaver_stats.Fmag.diff.mean(not_nan) + cadaver_stats.Fmag.diff.std(not_nan), fliplr(cadaver_stats.Fmag.diff.mean(not_nan) - cadaver_stats.Fmag.diff.std(not_nan))],...
     'k', 'FaceAlpha',alpha_std, 'EdgeColor','none');

% mark t-test significant points
scatter(cadaver_stats.Fmag.bins(cadaver_stats.Fmag.diff.h), cadaver_stats.Fmag.diff.mean(cadaver_stats.Fmag.diff.h), 40, 'm', 'o', 'LineWidth',1.7)

xlabel('Angular Insertion Depth (\circ)'); 
ylabel('\Delta ||F|| (mN)');

linkaxes(h_ax_t, 'x');
ylim([min(cadaver_stats.Fmag.diff.mean(not_nan) - cadaver_stats.Fmag.diff.std(not_nan) - 10), max(cadaver_stats.Fmag.diff.mean(not_nan) + cadaver_stats.Fmag.diff.std(not_nan) + 10)]) 
