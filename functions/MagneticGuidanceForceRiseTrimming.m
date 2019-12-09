function [ind_trimmed, X_cutoff, ind_X_cutoff] = MagneticGuidanceForceRiseTrimming(X, F, dF_thresh, X_span)
%% Determines where forces increase by a given magnitude within a specified span
%
%   X        => array of values corresponding to each value in F (e.g. linear depth)
%   F        => 1D array of forces
%   F_thresh => magnitude F must rise for cutoff
%   X_span   => span to check for F_rise (NOTE: same units as X)
%
%   ind_trimmed  => [size of X] logical indices corresponding to trimmed data -> X(ind_trimmed)
%   X_cutoff     => raw value of X at the cutoff point
%   ind_X_cutoff => index of X_cutoff -> X(ind_X_cutoff) = X_cutoff
%
%
%   Trevor Bruns
%   December 2019



%%
ind_trimmed = true(size(X));
X_cutoff = nan;
ind_X_cutoff = nan;

for ii=1:length(X)
    
    ind_next = find(X - X(ii) > X_span, 1);
    
    if isempty(ind_next)
        break;
    else
        if (F(ind_next) - F(ii)) > dF_thresh
            ind_X_cutoff = ii;
            X_cutoff = X(ind_X_cutoff);
            ind_trimmed(ind_X_cutoff:end) = false;
            break;
        end
    end
            
end

end