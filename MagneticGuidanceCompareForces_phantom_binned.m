%% Compare binned phantom forces

% Output is a binned plot and we will save the binned data for use with the
% colorbar plot

% Trevor Bruns and Katy Riojas
% Last Updated: December 2019


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equal Width Bins (non-normalized) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%
% Normalized Bins %
%%%%%%%%%%%%%%%%%%%

% before basal turn
for i_bin = 1:n_pre_bins
    
    % create vectors containing all force measurements within the corresponding bin
    nomag_normbin(i_bin).Fmags = [];
      mag_normbin(i_bin).Fmags = [];

    for i_trial = 1:length(data_robotic_phantom)
        % append forces from current trial
        nomag_normbin(i_bin).Fmags = [nomag_normbin(i_bin).Fmags; data_robotic_phantom(i_trial).nomag_normbin_Fmag.pre.vals{i_bin}];
          mag_normbin(i_bin).Fmags = [mag_normbin(i_bin).Fmags;   data_robotic_phantom(i_trial).mag_normbin_Fmag.pre.vals{i_bin}];
    end

    manual_normbin(i_bin).Fmags = [];
    for i_trial = 1:length(data_manual_phantom)
        % append forces from current trial
        manual_normbin(i_bin).Fmags = [manual_normbin(i_bin).Fmags; data_manual_phantom(i_trial).normbin_Fmag.pre.vals{i_bin}];
    end
    
    % compute mean for current bin
    phantom_stats.Fmag.normbin.mean.nomag(i_bin)  = mean(nomag_normbin(i_bin).Fmags);
    phantom_stats.Fmag.normbin.mean.mag(i_bin)    = mean(mag_normbin(i_bin).Fmags);
    phantom_stats.Fmag.normbin.mean.manual(i_bin) = mean(manual_normbin(i_bin).Fmags);

    % compute standard deviation
    phantom_stats.Fmag.normbin.std.nomag(i_bin)  = std(nomag_normbin(i_bin).Fmags);
    phantom_stats.Fmag.normbin.std.mag(i_bin)    = std(mag_normbin(i_bin).Fmags);
    phantom_stats.Fmag.normbin.std.manual(i_bin) = std(manual_normbin(i_bin).Fmags);

    % perform t-test to determine if mean force with magnet is less than without
    [phantom_stats.Fmag.normbin.diff.h(i_bin),  phantom_stats.Fmag.normbin.diff.p(i_bin), phantom_stats.Fmag.normbin.diff.ci(i_bin,:)] = ...
            ttest2( mag_normbin(i_bin).Fmags, nomag_normbin(i_bin).Fmags, 'Tail','left', 'Vartype','unequal', 'Alpha', 0.01);

    % compute mean difference between manual/robotic
    phantom_stats.Fmag.normbin.diff.mean(i_bin) = phantom_stats.Fmag.normbin.std.mag(i_bin) - phantom_stats.Fmag.normbin.std.nomag(i_bin);

    % compute RMS standard deviation of the difference
    phantom_stats.Fmag.normbin.diff.std(i_bin) = sqrt( phantom_stats.Fmag.normbin.std.nomag(i_bin)^2 + phantom_stats.Fmag.normbin.std.mag(i_bin)^2);

end

% after basal turn
for i_post_bin = 1:n_post_bins
    
    i_bin = n_pre_bins + i_post_bin;

    % create vectors containing all force measurements within the corresponding bin
    nomag_normbin(i_bin).Fmags = [];
      mag_normbin(i_bin).Fmags = [];

    for i_trial = 1:length(data_robotic_phantom)
        % append forces from current trial
        nomag_normbin(i_bin).Fmags = [nomag_normbin(i_bin).Fmags; data_robotic_phantom(i_trial).nomag_normbin_Fmag.post.vals{i_post_bin}];
          mag_normbin(i_bin).Fmags = [mag_normbin(i_bin).Fmags;   data_robotic_phantom(i_trial).mag_normbin_Fmag.post.vals{i_post_bin}];
    end

    manual_normbin(i_bin).Fmags = [];
    for i_trial = 1:length(data_manual_phantom)
        % append forces from current trial
        manual_normbin(i_bin).Fmags = [manual_normbin(i_bin).Fmags; data_manual_phantom(i_trial).normbin_Fmag.post.vals{i_post_bin}];
    end
    
    % compute mean for current bin
    phantom_stats.Fmag.normbin.mean.nomag(i_bin)  = mean(nomag_normbin(i_bin).Fmags);
    phantom_stats.Fmag.normbin.mean.mag(i_bin)    = mean(mag_normbin(i_bin).Fmags);
    phantom_stats.Fmag.normbin.mean.manual(i_bin) = mean(manual_normbin(i_bin).Fmags);

    % compute standard deviation
    phantom_stats.Fmag.normbin.std.nomag(i_bin)  = std(nomag_normbin(i_bin).Fmags);
    phantom_stats.Fmag.normbin.std.mag(i_bin)    = std(mag_normbin(i_bin).Fmags);
    phantom_stats.Fmag.normbin.std.manual(i_bin) = std(manual_normbin(i_bin).Fmags);

    % perform t-test to determine if mean force with magnet is less than without
    [phantom_stats.Fmag.normbin.diff.h(i_bin),  phantom_stats.Fmag.normbin.diff.p(i_bin), phantom_stats.Fmag.normbin.diff.ci(i_bin,:)] = ...
            ttest2( mag_normbin(i_bin).Fmags, nomag_normbin(i_bin).Fmags, 'Tail','left', 'Vartype','unequal', 'Alpha', 0.01);

    % compute mean difference between manual/robotic
    phantom_stats.Fmag.normbin.diff.mean(i_bin) = phantom_stats.Fmag.normbin.mean.mag(i_bin) - phantom_stats.Fmag.normbin.mean.nomag(i_bin);

    % compute RMS standard deviation of the difference
    phantom_stats.Fmag.normbin.diff.std(i_bin) = sqrt( phantom_stats.Fmag.normbin.std.nomag(i_bin)^2 + phantom_stats.Fmag.normbin.std.mag(i_bin)^2);

end

% remove NaNs and convert to logicals
phantom_stats.Fmag.normbin.diff.h(isnan(phantom_stats.Fmag.normbin.diff.h)) = 0;
phantom_stats.Fmag.normbin.diff.h = logical(phantom_stats.Fmag.normbin.diff.h);



%% Fmag vs Normalized Insertion Depth (all trials)
figure(18); clf(18);
hold on; grid on;
title('Phantom: ||F|| vs Normalized Insertion Depth');

line_width = 1;

for i_trial = 1:length(data_manual_phantom)
    % manual
    plot([data_manual_phantom(i_trial).normbin_Fmag.pre.mean; data_manual_phantom(i_trial).normbin_Fmag.post.mean],...
        'Color', 'r', 'LineWidth', line_width);
    % unguided
    plot([data_robotic_phantom(i_trial).nomag_normbin_Fmag.pre.mean; data_robotic_phantom(i_trial).nomag_normbin_Fmag.post.mean],...
        'Color', 'b', 'LineWidth', line_width);
    % guided
    plot([data_robotic_phantom(i_trial).mag_normbin_Fmag.pre.mean; data_robotic_phantom(i_trial).mag_normbin_Fmag.post.mean],...
        'Color', 'g', 'LineWidth', line_width);
end

% mark basal turn point
xline(n_pre_bins, '--k', 'Basal Turn', 'LineWidth',2, 'LabelHorizontalAlignment','center', 'LabelVerticalAlignment','middle');

xlabel('Normalized Angular Insertion Depth');
ylabel('||F|| (mN)');
set(gca,'xticklabel',{[]})

legend('Manual','Robotic','Robotic & Magnetic Steering', 'Location','nw');



%% Fmag vs Normalized Insertion Depth (averages)

if exist('hf_normbin','var')
    if isvalid(hf_normbin)
        close(hf_normbin)
    end
end

hf_normbin = figure;
hf_normbin.WindowState = "maximized";

line_width = 2;
alpha_std = 0.15;

% h_ax(1) = subplot_er(2,1,1);
h_ax(1) = subplot(2,1,1);
grid on; hold on; %xlim([0 400]);
title(sprintf('Phantom: ||F|| vs Normalized Insertion Depth (pre/post basal turn bins = [%i, %i]', n_pre_bins, n_post_bins));


% plot means
plot(phantom_stats.Fmag.normbin.mean.manual, 'Color', 'r','LineWidth',line_width);
plot(phantom_stats.Fmag.normbin.mean.nomag,  'Color', 'b','LineWidth',line_width);
plot(phantom_stats.Fmag.normbin.mean.mag,    'Color', 'g','LineWidth',line_width);

% plot standard deviations as shaded regions around means
n_normbins = n_pre_bins + n_post_bins;
X = 1:n_normbins;

not_nan = ~isnan(phantom_stats.Fmag.normbin.mean.manual);
fill([X(not_nan), fliplr(X(not_nan))],...
     [phantom_stats.Fmag.normbin.mean.manual(not_nan) + phantom_stats.Fmag.normbin.std.manual(not_nan), fliplr(phantom_stats.Fmag.normbin.mean.manual(not_nan) - phantom_stats.Fmag.normbin.std.manual(not_nan))],...
     'r', 'FaceAlpha',alpha_std, 'EdgeColor','none');

not_nan = ~isnan(phantom_stats.Fmag.normbin.mean.nomag);
fill([X(not_nan), fliplr(X(not_nan))],...
     [phantom_stats.Fmag.normbin.mean.nomag(not_nan) + phantom_stats.Fmag.normbin.std.nomag(not_nan), fliplr(phantom_stats.Fmag.normbin.mean.nomag(not_nan) - phantom_stats.Fmag.normbin.std.nomag(not_nan))],...
     'b', 'FaceAlpha',alpha_std, 'EdgeColor','none');

not_nan = ~isnan(phantom_stats.Fmag.normbin.mean.mag);
fill([X(not_nan), fliplr(X(not_nan))],...
     [phantom_stats.Fmag.normbin.mean.mag(not_nan) + phantom_stats.Fmag.normbin.std.mag(not_nan), fliplr(phantom_stats.Fmag.normbin.mean.mag(not_nan) - phantom_stats.Fmag.normbin.std.mag(not_nan))],...
     'g', 'FaceAlpha',alpha_std, 'EdgeColor','none');


ylabel('||F|| (mN)');
legend('Manual','Robotic','Robotic & Magnetic Steering', 'Location','nw');
ylim([-10, h_ax(1).YLim(2)])

% mark basal turn point
xline(n_pre_bins, '--k', 'Basal Turn', 'LineWidth',2, 'LabelHorizontalAlignment','center', 'LabelVerticalAlignment','middle');
set(gca,'xticklabel',{[]})



% Plot delta F between guided/unguided robotic insertion means
h_ax(2) = subplot(2,1,2);
% h_ax(2) = subplot_er(2,1,2);
grid on; hold on;

% delta F
plot(phantom_stats.Fmag.normbin.diff.mean, 'Color', 'k','LineWidth',line_width);

% standard deviation
not_nan = ~isnan(phantom_stats.Fmag.normbin.diff.mean);
fill([X(not_nan), fliplr(X(not_nan))],...
     [phantom_stats.Fmag.normbin.diff.mean(not_nan) + phantom_stats.Fmag.normbin.diff.std(not_nan), fliplr(phantom_stats.Fmag.normbin.diff.mean(not_nan) - phantom_stats.Fmag.normbin.diff.std(not_nan))],...
     'k', 'FaceAlpha',alpha_std, 'EdgeColor','none');

% mark t-test significant points
scatter(X(phantom_stats.Fmag.normbin.diff.h), phantom_stats.Fmag.normbin.diff.mean(phantom_stats.Fmag.normbin.diff.h), 40, 'm', 'o', 'LineWidth',1.7)

% mark basal turn point
xline(n_pre_bins, '--k', 'Basal Turn', 'LineWidth',2, 'LabelHorizontalAlignment','center', 'LabelVerticalAlignment','middle');

xlabel('Normalized Angular Insertion Depth');
ylabel('\Delta ||F|| (mN)');
set(gca,'xticklabel',{[]})

linkaxes(h_ax, 'x');


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
grid on; hold on; %xlim([0 400]);
title(strcat(sprintf('Mean Forces in Phantom (Bin Size = %i', phantom_stats.Fmag.bins(2)-phantom_stats.Fmag.bins(1)), '\circ)'))

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
    scatter(phantom_stats.Fmag.bins(last_ind), phantom_stats.Fmag.mean.manual(last_ind), 120, 'r', 'd', 'filled');

    last_ind = data_robotic_phantom(i_trial).nomag_mea_binned.ind(end);
    scatter(phantom_stats.Fmag.bins(last_ind), phantom_stats.Fmag.mean.nomag(last_ind),  120, 'b', 'd', 'filled');

    last_ind = data_robotic_phantom(i_trial).mag_binned.ind(end-1); % TODO: fix NaN in interp_angdepth
    scatter(phantom_stats.Fmag.bins(last_ind), phantom_stats.Fmag.mean.mag(last_ind),    120, 'g', 'd', 'filled');
end

ylabel('||F|| (mN)');
legend('Manual','Robotic','Robotic & Magnetic Steering', 'Location','nw');
ylim([-10, h_ax(1).YLim(2)])
% ylim([0 120])

% Plot delta F between guided/unguided robotic insertion means
h_ax(2) = subplot_er(2,1,2);
grid on; hold on;

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
scatter(phantom_stats.Fmag.bins(phantom_stats.Fmag.diff.h), phantom_stats.Fmag.diff.mean(phantom_stats.Fmag.diff.h), 40, 'm', 'o', 'LineWidth',1.7)

xlabel('Angular Insertion Depth (\circ)'); 
ylabel('\Delta ||F|| (mN)');

linkaxes(h_ax, 'x');
xlim([0, phantom_stats.Fmag.bins(end)+10])
% xlim([0,400])
ylim([min(phantom_stats.Fmag.diff.mean(not_nan) - phantom_stats.Fmag.diff.std(not_nan) - 10), max(phantom_stats.Fmag.diff.mean(not_nan) + phantom_stats.Fmag.diff.std(not_nan) + 10)]) 


%% Plot p-values from t-test
% figure(20); clf(20);
% bar(phantom_stats.Fmag.bins, phantom_stats.Fmag.diff.p)


%% Update Saved Phantom Data if it is called for
% if update_saved_phantom_structs
%    save('data\phantom\data_manual_phantom.mat','data_manual_phantom');
%    save('data\phantom\data_robotic_phantom.mat','data_robotic_phantom');
% end