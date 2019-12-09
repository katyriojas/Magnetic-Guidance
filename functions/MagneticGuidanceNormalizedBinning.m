function norm_bin = MagneticGuidanceNormalizedBinning(X, Y, threshold, n_pre_bins, n_post_bins)
%% Separates X into two partitions and bins the corresponding Y values into the specified number of pre and post bins
%
% Example: 
%   X = interp_angdepth
%   Y = Fmag
%   threshold   = basal_turn_angle
%   n_pre_bins  = 40
%   n_post_bins = 120


%% Pre-Threshold

% indices in X less than threshold
X_pre = find(X<=threshold);

% create vector of the bin edges
edges = linspace(0, threshold, n_pre_bins+1);

% create vector of the bin centers
bin_width = edges(2) - edges(1);
bins = edges(2:end) - bin_width/2;

% sort X(pre_ind) into bins
[~, ~, bin_ind] = histcounts(X(X_pre), edges);

% sort corresponding Y values into bins and compute mean/std
vals  = cell(n_pre_bins,1);
means  = nan(n_pre_bins, 1);
stdev = nan(n_pre_bins, 1);

for i_bin = 1:n_pre_bins

    vals{i_bin}  = Y( X_pre(bin_ind == i_bin) );

    if ~isempty(vals{i_bin})
        means(i_bin)  = mean(vals{i_bin});
        stdev(i_bin) = std(vals{i_bin});
    end
end

% save into norm_bin struct
norm_bin.pre.bins  = bins;
norm_bin.pre.edges = edges;
norm_bin.pre.vals  = vals;
norm_bin.pre.mean  = means;
norm_bin.pre.stdev = stdev;


%% Post-Threshold

% indices in X greater than threshold
X_post = find(X>threshold);

% create vector of the bin edges
edges = linspace(threshold, max(X(X_post)), n_post_bins+1);

% create vector of the bin centers
bin_width = edges(2) - edges(1);
bins = edges(2:end) - bin_width/2;

% sort X(post_ind) into bins
[~, ~, bin_ind] = histcounts(X(X_post), edges);

% sort corresponding Y values into bins and compute mean/std
vals  = cell(n_post_bins,1);
means  = nan(n_post_bins, 1);
stdev = nan(n_post_bins, 1);
for i_bin = 1:n_post_bins

    vals{i_bin}  = Y( X_post(bin_ind == i_bin) );

    if ~isempty(vals{i_bin})
        means(i_bin)  = mean(vals{i_bin});
        stdev(i_bin) = std(vals{i_bin});
    end
end

% save into norm_bin struct
norm_bin.post.bins  = bins;
norm_bin.post.edges = edges;
norm_bin.post.vals  = vals;
norm_bin.post.mean  = means;
norm_bin.post.stdev = stdev;

end