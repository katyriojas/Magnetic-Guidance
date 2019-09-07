function [scala_tympani_path, theta] = scalaTympaniMedialAxis(varargin)
%   Outputs the medial axis for a typical human scala tympani.
%   
%   Equations and constants from:
%   A Scalable Model for Human Scala-Tympani Phantoms (Clark 2011)
%   &
%   Scala-Tympani Phantom With Cochleostomy and Round-Window Openings for 
%   Cochlear-Implant Insertion Experiments (Leon 2014)
%   
%
%   INPUTS (optional)
%   theta_step => [degrees] angular spacing to use between points (default=0.1)
%   plot_flag  => [bool] if true, will create plots useful for debugging (default = false)
%
%   OUTPUTS
%   scala_tympani => [mm] 3xN Cartesian coordiates in standard clinical frame
%   theta => [deg] corresponding angle to each point in 
%
%   EXAMPLES:
%   scala_tympani = scalaTympaniMedialAxis();
%   use default spacing of theta=0.1 deg; no plots
%   
%   scala_tympani = scalaTympaniMedialAxis(5000, true);
%   scala_tympani will be [3x5000]; will create plots
%
%
%   Trevor Bruns
%   July 2019

%% parse inputs

plot_flag = false;
theta_step = 0.1;

if nargin > 0
    for ii=1:length(varargin)
        if isa(varargin{ii},'double')
            theta_step = varargin{ii};
        elseif isa(varargin{ii},'logical')
            plot_flag = varargin{ii};
        end
    end
else
    error('too many inputs')        
end

%% constants
A = 3.762;     % [mm]
B = 0.001317;  % [mm]
C = 7.967;     % [mm]
D = 0.1287;    % [mm]
E = 0.003056;  % [mm]
theta0 = 5.0;  % [deg]
theta1 = 10.3; % [deg]

%% base curve equations (in cylindrical coordinates)
theta_start = 6.8; % [deg]
theta_end = 910.3; % [deg]

theta = theta_start:theta_step:theta_end; % [deg] angular range of ST curve
R = zeros(1,length(theta)); % [mm] radial distance from spiral center
Z = zeros(1,length(theta)); % [mm] height

for ii=1:length(theta)
   if theta(ii) < 100
       R(ii) = C*(1-D*log(theta(ii)-theta0));
   else
       R(ii) = A*exp(-B*(0.0002*theta(ii)^2 + 0.98*theta(ii)));
   end
   
   Z(ii) = E*(theta(ii)-theta1);
end

%% smooth to remove sharp transitions

% first want to remove transition at 100 degrees
smooth_span = 0.020; % [fraction] 1.0 = all points
R_smooth_full = smooth(theta, R, smooth_span, 'loess')';

% smooth distorts values at the beginning -> want to keep original values
k = find(theta>80,1); % index where theta = 80 degrees
R_smooth = R_smooth_full;
R_smooth(1:k) = R(1:k);

%% convert to Cartesian coordinates
scala_tympani_orig = [R.*cosd(theta); R.*sind(theta); Z];
scala_tympani_smooth = [R_smooth.*cosd(theta); R_smooth.*sind(theta); Z];

%% compute medial shift adjustments (as shown in Fig 8)

% constants (values from Table 3)
% Note: added a repeated value at end for 910.3 deg to prevent NaN when interpolating
theta_s = [ 6.8, 9.5, 15.1, 19.0, 23.5, 34.3, 46.7, 60.4, 75.3, 91.3, 108.4, 126.1 ...
          , 144.2, 162.8, 181.9, 201.5, 221.6, 242.4, 263.7, 285.8, 308.5, 332.0 ...
          , 356.4, 381.6, 407.8, 435.1, 463.6, 493.4, 524.6, 557.3, 591.9, 628.4 ...
          , 667.2, 708.6, 753.0, 800.9, 853.0, 910.1, 910.3];

d = [0, 1, 2, 2.5, 3:36, 36.02]; 
      
phi_s = [ 1.57, 1.57, 1.43, 1.17, 0.97, 0.56, 0.00, 0.02, 0.04, 0.06:0.01:0.30 ...
        , 0.25, 0.20, 0.15, 0.10, 0.10 ]; % [rad]
    
w_s = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.3, 0.28, 0.26, 0.24, 0.28 ...
     , 0.28, 0.32, 0.35, 0.40, 0.50, 0.55, 0.62, 0.62, 0.60, 0.60 ...
     , 0.65, 0.66, 0.67, 0.68, 0.65, 0.65, 0.65, 0.57, 0.49, 0.41 ...
     , 0.33, 0.25, 0.17, 0.09, 0.01, 0.01, 0.01 ]; % [mm]
 
% compute cumulative length along scala tympani (use orig to ensure we match with published data)
d_orig = [0, cumsum(sqrt(sum(diff(scala_tympani_orig')'.^2,1)))];

% interpolate to match values with d, then smooth slightly
phi_s_interp = interp1(d, phi_s, d_orig);
phi_s_smooth = smooth(d_orig, phi_s_interp, 0.001, 'lowess')';

w_s_interp = interp1(d, w_s, d_orig);
w_s_smooth = smooth(d_orig, w_s_interp, 0.001, 'lowess')';


% % interpolate to match values with theta
% phi_s_interp = interp1(theta_s, phi_s, theta);
% w_s_interp = interp1(theta_s, w_s, theta);

% compute adjustments in local frame [x',z']
adj_x_prime = -w_s_smooth.*cos(phi_s_smooth);
adj_z_prime = -w_s_smooth.*sin(phi_s_smooth);


%% transform adjustments into global frame

% compute alpha => angle of normal vector w.r.t. x-axis
ST_gradient = gradient(scala_tympani_orig);
ST_gradient_smooth = gradient(scala_tympani_smooth);
alpha = bsxfun(@atan2d, ST_gradient(2,:), ST_gradient(1,:)) - 90;
alpha_smooth = bsxfun(@atan2d, ST_gradient_smooth(2,:), ST_gradient_smooth(1,:)) - 90;

% x' unit vector at each point in cartesian coordinates
x_prime = [cosd(alpha_smooth); sind(alpha_smooth)];

%% apply adjustments
scala_tympani_adj = [scala_tympani_smooth(1,:) + adj_x_prime.*x_prime(1,:); ...
                 scala_tympani_smooth(2,:) + adj_x_prime.*x_prime(2,:); ...
                 scala_tympani_smooth(3,:) + adj_z_prime];

%% choose path to output
scala_tympani_path = scala_tympani_smooth;

%% plot
if plot_flag
    figure(1)
    clf(1)
    plot(theta,alpha,'k--'); 
    hold on; 
    plot(theta,alpha_smooth,'-.b'); 
    legend('alpha','alpha\_smooth')

    figure(2)
    clf(2)
    plot(theta,R,'k--'); 
    hold on; 
    plot(theta,R_smooth,'-.b'); 
    legend('R','R\_smooth')

    figure(3)
    clf(3)
    drawPolyline3d(scala_tympani_orig', ':k');
    hold on;
    drawPolyline3d(scala_tympani_smooth', '--r');
    drawPolyline3d(scala_tympani', 'b');
    legend('base','smooth','smooth\_adj')
    axis equal
end
end