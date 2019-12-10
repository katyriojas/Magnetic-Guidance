%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Performs all binning and t-test computations for RA-L data            %
%                                                                         %
%   Must first run LoadRALData_Robotic_Phantom.m & LoadRALData_Manual.m   %
%                                                                         %
%   Trevor Bruns & Katy Riojas                                            %
%   December 2019                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

degspan = 3; % [deg] span/width of each bin

%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   %%%
%%%  Robotic Phantom  %%%
%%%                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Binning - Angular Depth (untrimmed)

% compute the largest angular insertion depth to bound our bins
max_ang = 0;
for ii = 1:length(data_robotic_phantom)
    max_ang = max([max_ang,...
                   max(data_robotic_phantom(ii).nomag_mea_interp_angdepth),...
                   max(data_robotic_phantom(ii).mag_interp_angdepth)]);
end
num_bins = ceil(max_ang/degspan);

% create vector of the bin edges
bin_edges = 0:degspan:num_bins*degspan;

% create vector of the bin centers
bins = bin_edges(2:end) - degspan/2;

% sort into bins
for ii = 1:length(data_robotic_phantom)

    % determine the bin for each interp_angdepth point
    [~,~,data_robotic_phantom(ii).nomag_mea_binned.ind] = histcounts(data_robotic_phantom(ii).nomag_mea_interp_angdepth, bin_edges);
    [~,~,data_robotic_phantom(ii).mag_binned.ind]       = histcounts(data_robotic_phantom(ii).mag_interp_angdepth, bin_edges);

    % compute mean of all the force measurements within each bin (NaN for bins with no measurements)
    data_robotic_phantom(ii).nomag_mea_binned.Fmean = zeros(size(bins));
    data_robotic_phantom(ii).mag_binned.Fmean = zeros(size(bins));
    for jj = 1:length(bins)
        data_robotic_phantom(ii).nomag_mea_binned.Fmean(jj) = ...
            mean( data_robotic_phantom(ii).nomag_mea.Fmag( data_robotic_phantom(ii).nomag_mea_binned.ind == jj ) );

        data_robotic_phantom(ii).mag_binned.Fmean(jj) = ...
            mean( data_robotic_phantom(ii).mag.Fmag( data_robotic_phantom(ii).mag_binned.ind == jj ) );
    end
    

    % also save the actual bins and edges
    data_robotic_phantom(ii).nomag_mea_binned.bins  = bins;
    data_robotic_phantom(ii).nomag_mea_binned.edges = bin_edges;
    data_robotic_phantom(ii).mag_binned.bins  = bins;
    data_robotic_phantom(ii).mag_binned.edges = bin_edges;
end



%% Binning - Angular Depth (trimmed)

% compute the largest angular insertion depth to bound our bins
max_ang = 0;
for ii = 1:length(data_robotic_phantom_trim)
    max_ang = max([max_ang,...
                   max(data_robotic_phantom_trim(ii).nomag_mea_interp_angdepth),...
                   max(data_robotic_phantom_trim(ii).mag_interp_angdepth)]);
end
num_bins = ceil(max_ang/degspan);

% create vector of the bin edges
bin_edges = 0:degspan:num_bins*degspan;

% create vector of the bin centers
bins = bin_edges(2:end) - degspan/2;

% sort into bins
for ii = 1:length(data_robotic_phantom_trim)

    % determine the bin for each interp_angdepth point
    [~,~,data_robotic_phantom_trim(ii).nomag_mea_binned.ind] = histcounts(data_robotic_phantom_trim(ii).nomag_mea_interp_angdepth, bin_edges);
    [~,~,data_robotic_phantom_trim(ii).mag_binned.ind]       = histcounts(data_robotic_phantom_trim(ii).mag_interp_angdepth, bin_edges);

    % compute mean of all the force measurements within each bin (NaN for bins with no measurements)
    data_robotic_phantom_trim(ii).nomag_mea_binned.Fmean = zeros(size(bins));
    data_robotic_phantom_trim(ii).mag_binned.Fmean = zeros(size(bins));
    for jj = 1:length(bins)
        data_robotic_phantom_trim(ii).nomag_mea_binned.Fmean(jj) = ...
            mean( data_robotic_phantom_trim(ii).nomag_mea.Fmag( data_robotic_phantom_trim(ii).nomag_mea_binned.ind == jj ) );

        data_robotic_phantom_trim(ii).mag_binned.Fmean(jj) = ...
            mean( data_robotic_phantom_trim(ii).mag.Fmag( data_robotic_phantom_trim(ii).mag_binned.ind == jj ) );
    end

    % also save the actual bins and edges
    data_robotic_phantom_trim(ii).nomag_mea_binned.bins  = bins;
    data_robotic_phantom_trim(ii).nomag_mea_binned.edges = bin_edges;
    data_robotic_phantom_trim(ii).mag_binned.bins  = bins;
    data_robotic_phantom_trim(ii).mag_binned.edges = bin_edges;
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  %%%
%%%  Manual Phantom  %%%
%%%                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%


%% Binning

% sort into bins
for ii = 1:length(data_manual_phantom)

    % determine the bin for each interp_angdepth point
    [~,~,data_manual_phantom(ii).binned.ind] = histcounts(data_manual_phantom(ii).interp_angdepth, data_robotic_phantom(ii).mag_binned.edges);

    % compute mean of all the force measurements within each bin (NaN for bins with no measurements)
    data_manual_phantom(ii).binned.Fmean = zeros(size(bins));
    for jj = 1:length(bins)
        data_manual_phantom(ii).binned.Fmean(jj) = ...
            mean( data_manual_phantom(ii).Fmag_trimmed( data_manual_phantom(ii).binned.ind == jj ) );
    end
    
    % also save the actual bins and edges used
    data_manual_phantom(ii).binned.bins  = bins;
    data_manual_phantom(ii).binned.edges = bin_edges;
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        %%%
%%%  Statistics (trimmed)  %%%
%%%                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Combine binned force measurements for each set of trials and compute statistics

clear phantom_stats_trim;

phantom_stats_trim.Fmag.bins = data_robotic_phantom_trim(1).mag_binned.bins;

for i_bin = 1:length(phantom_stats_trim.Fmag.bins)
    
    % create vectors containing all force measurements within the corresponding bin
    nomag_mea_binned(i_bin).Fmags = [];
          mag_binned(i_bin).Fmags = [];

    for i_trial = 1:length(data_robotic_phantom_trim)
        % append forces from current trial
        nomag_mea_binned(i_bin).Fmags = [nomag_mea_binned(i_bin).Fmags; ...
                                      data_robotic_phantom_trim(i_trial).nomag_mea.Fmag(data_robotic_phantom_trim(i_trial).nomag_mea_binned.ind == i_bin)];

        mag_binned(i_bin).Fmags = [mag_binned(i_bin).Fmags; ...
                                data_robotic_phantom_trim(i_trial).mag.Fmag(data_robotic_phantom_trim(i_trial).mag_binned.ind == i_bin)];
    end

    manual_binned(i_bin).Fmags = [];
    for i_trial = 1:length(data_manual_phantom)
        % append forces from current trial
        manual_binned(i_bin).Fmags = [manual_binned(i_bin).Fmags; ...
                                      data_manual_phantom(i_trial).Fmag_trimmed(data_manual_phantom(i_trial).binned.ind == i_bin)];
    end
    
    % compute mean for current bin
    phantom_stats_trim.Fmag.mean.nomag(i_bin)  = mean(nomag_mea_binned(i_bin).Fmags);
    phantom_stats_trim.Fmag.mean.mag(i_bin)    = mean(mag_binned(i_bin).Fmags);
    phantom_stats_trim.Fmag.mean.manual(i_bin) = mean(manual_binned(i_bin).Fmags);

    % compute standard deviation
    phantom_stats_trim.Fmag.std.nomag(i_bin)  = std(nomag_mea_binned(i_bin).Fmags);
    phantom_stats_trim.Fmag.std.mag(i_bin)    = std(mag_binned(i_bin).Fmags);
    phantom_stats_trim.Fmag.std.manual(i_bin) = std(manual_binned(i_bin).Fmags);

    % perform t-test to determine if mean force with magnet is less than without
    [phantom_stats_trim.Fmag.diff.h(i_bin),  phantom_stats_trim.Fmag.diff.p(i_bin), phantom_stats_trim.Fmag.diff.ci(i_bin,:)] = ...
            ttest2( mag_binned(i_bin).Fmags, nomag_mea_binned(i_bin).Fmags, 'Tail','left', 'Vartype','unequal', 'Alpha', 0.05);

    % compute mean difference between manual/robotic
    phantom_stats_trim.Fmag.diff.mean(i_bin) = phantom_stats_trim.Fmag.mean.mag(i_bin) - phantom_stats_trim.Fmag.mean.nomag(i_bin);

    % compute standard error of the difference between means
    phantom_stats_trim.Fmag.diff.std(i_bin) = sqrt( (phantom_stats_trim.Fmag.std.nomag(i_bin)^2) / length(nomag_mea_binned(i_bin).Fmags) ...
                                                    + (phantom_stats_trim.Fmag.std.mag(i_bin)^2) / length(mag_binned(i_bin).Fmags) );
end

% remove NaNs and convert to logicals
phantom_stats_trim.Fmag.diff.h(isnan(phantom_stats_trim.Fmag.diff.h)) = 0;
phantom_stats_trim.Fmag.diff.h = logical(phantom_stats_trim.Fmag.diff.h);

% peak forces (binned)
for i_trial = 1:length(data_robotic_phantom_trim)
    data_robotic_phantom_trim(i_trial).mag_binned.Fmax       = max(data_robotic_phantom_trim(i_trial).mag_binned.Fmean);
    data_robotic_phantom_trim(i_trial).nomag_mea_binned.Fmax = max(data_robotic_phantom_trim(i_trial).nomag_mea_binned.Fmean);
    data_manual_phantom(i_trial).binned.Fmax                 = max(data_manual_phantom(i_trial).binned.Fmean);
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          %%%
%%%  Statistics (untrimmed)  %%%
%%%                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
            ttest2( mag_binned(i_bin).Fmags, nomag_mea_binned(i_bin).Fmags, 'Tail','left', 'Vartype','unequal', 'Alpha', 0.05);

    % compute mean difference between manual/robotic
    phantom_stats.Fmag.diff.mean(i_bin) = phantom_stats.Fmag.mean.mag(i_bin) - phantom_stats.Fmag.mean.nomag(i_bin);

    % compute standard error of the difference between means
    phantom_stats.Fmag.diff.std(i_bin) = sqrt( (phantom_stats.Fmag.std.nomag(i_bin)^2) / length(nomag_mea_binned(i_bin).Fmags) ...
                                               + (phantom_stats.Fmag.std.mag(i_bin)^2) / length(mag_binned(i_bin).Fmags) );

end

% remove NaNs and convert to logicals
phantom_stats.Fmag.diff.h(isnan(phantom_stats.Fmag.diff.h)) = 0;
phantom_stats.Fmag.diff.h = logical(phantom_stats.Fmag.diff.h);

% peak forces (binned)
for i_trial = 1:length(data_robotic_phantom_trim)
    data_robotic_phantom(i_trial).mag_binned.Fmax = max(data_robotic_phantom(i_trial).mag_binned.Fmean);
    data_robotic_phantom(i_trial).nomag_mea_binned.Fmax = max(data_robotic_phantom(i_trial).nomag_mea_binned.Fmean);
end
