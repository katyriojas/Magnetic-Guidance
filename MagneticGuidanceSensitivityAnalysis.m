%% Evaluates sensitivity of magnetic field to errors in Omnimagnet position and magnitude
% 
% Uses Jake Abbott's code, but evaluates at all points along ST path used for RA-L experiments
%
% Trevor Bruns/Jake Abbott
% December 2019
%

%% Magnetic field

%%
%%%%%%%%%%%
% Phantom %
%%%%%%%%%%%

% Import points along ST path (p_stpath) and their associated field vectors (b_stpath)
load('data/phantom/phantom_RAL_sensitivityparams.mat');

% trim away points where magnet is turned off
i_on  = find( abs(b_stpath(1,:)) > 0, 1);
i_off = i_on-2 + find( abs(b_stpath(1,i_on:end)) < 1e-6, 1);
b_nom = b_stpath(:, i_on:i_off) * 1e-3; % also convert from [mT] to [T]
p_nom = p_stpath(:, i_on:i_off) * 1e-3; % [mm] to [m]

% Compute max errors at each path point
for ii = 1:length(p_nom)
    [pos_error_mag(ii), pos_error_ang(ii), ang_error_mag(ii), ang_error_ang(ii)] = OmnimagnetSensitivity(p_nom(:,ii), b_nom(:,ii));
end

phantom_max_pos_error_mag = max(pos_error_mag)
phantom_pos_error_ang     = max(pos_error_ang)
phantom_ang_error_mag     = max(ang_error_mag)
phantom_ang_error_ang     = max(ang_error_ang)

%%
%%%%%%%%%%%
% Cadaver %
%%%%%%%%%%%

% Import points along ST path (p_stpath) and their associated field vectors (b_stpath)
load('data/cadaver/cadaver_RAL_sensitivityparams.mat');

% trim away points where magnet is turned off
i_on  = find( abs(b_stpath(1,:)) > 0, 1);
i_off = i_on-2 + find( abs(b_stpath(1,i_on:end)) < 1e-6, 1);
b_nom = b_stpath(:, i_on:i_off) * 1e-3; % also convert from [mT] to [T]
p_nom = p_stpath(:, i_on:i_off) * 1e-3; % [mm] to [m]

% Compute max errors at each path point
for ii = 1:length(p_nom)
    [pos_error_mag(ii), pos_error_ang(ii), ang_error_mag(ii), ang_error_ang(ii)] = OmnimagnetSensitivity(p_nom(:,ii), b_nom(:,ii));
end

cadaver_max_pos_error_mag = max(pos_error_mag)
cadaver_pos_error_ang     = max(pos_error_ang)
cadaver_ang_error_mag     = max(ang_error_mag)
cadaver_ang_error_ang     = max(ang_error_ang)