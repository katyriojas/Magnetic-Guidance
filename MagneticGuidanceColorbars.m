% Plot colorbar plot for all phantom data
% Output is a figure with every individual trial (manual, UG, G) plotted as
% a colorbar rectangle with the shaded region indicated Fmagsmooth

% Trevor Bruns and Katy Riojas
% Last Updated: 11/19/19

% TODO: Change so we don't need to run the bin script before running the
% colorbar script

%%
regenerate_manual_data = false;
regenerate_phantom_data = false;


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

%% Trim to 10 degrees from final angular depth for consistent comparison
N = size(data_robotic_phantom,2); % assumes equal # of manual, unguided, and guided trials

for ii = 1:N
    
    trim_ind_manual = find(data_manual_phantom(ii).binned.angle > (data_manual_phantom(ii).binned.angle(end)-10), 1);
    trim_ind_nomag  = find(data_robotic_phantom(ii).nomag_mea_binned.angle > (data_robotic_phantom(ii).nomag_mea_binned.angle(end)-10), 1);
    trim_ind_mag    = find(data_robotic_phantom(ii).mag_binned.angle > (data_robotic_phantom(ii).mag_binned.angle(end)-10), 1);
    
    phantom(ii).manual.angdepth = data_manual_phantom(ii).binned.angle(1:trim_ind_manual);
    phantom(ii).manual.Fmag_smooth = data_manual_phantom(ii).binned.Fmag(1:trim_ind_manual);

    phantom(ii).nomag.angdepth = data_robotic_phantom(ii).nomag_mea_binned.angle(1:trim_ind_nomag);
    phantom(ii).nomag.Fmag_smooth = data_robotic_phantom(ii).nomag_mea_binned.Fmag(1:trim_ind_nomag);

    phantom(ii).mag.angdepth = data_robotic_phantom(ii).mag_binned.angle(1:trim_ind_mag);
    phantom(ii).mag.Fmag_smooth = data_robotic_phantom(ii).mag_binned.Fmag(1:trim_ind_mag);
    
end


%% Create base profile (circle/square)
base_type = 'square';
width = 50;

if strcmp(base_type, 'square') 
    q = linspace(0,2*pi,5) + pi/4;
    base = (width/2)*[cos(q); sin(q)];   % Base curve is a square
    spacing = .7*width;
    cases_gap = 0.15*width;
elseif strcmp(base_type, 'circle')
    % Base curve is a circle
    cir_res = 40;   % circle resolution (# pts)
    q = linspace(0, 2*pi, cir_res);
    base = (width/2)*[cos(q); sin(q)];
    spacing = 1.02*width;
    cases_gap = 0.2*width;
else
    error('base_type must be ''circle'' or ''square''')
end

y_locations = linspace(3*spacing*N, 0, 3*N) - 1.5*spacing*N;
y_locations = y_locations + [cases_gap*ones(1,N), zeros(1,N), -cases_gap*ones(1,N)];
  
%% Create paths

for ii=1:N
    path(ii).manual = [phantom(ii).manual.angdepth,...
                       y_locations(ii)*ones(size(phantom(ii).manual.angdepth)),...
                       zeros(size(phantom(ii).manual.angdepth))]';

    path(ii).nomag  = [phantom(ii).nomag.angdepth,...
                       y_locations(ii+N)*ones(size(phantom(ii).nomag.angdepth)),...
                       zeros(size(phantom(ii).nomag.angdepth))]';

    path(ii).mag    = [phantom(ii).mag.angdepth,...
                      y_locations(ii+N*2)*ones(size(phantom(ii).mag.angdepth)),...
                      zeros(size(phantom(ii).mag.angdepth))]';
end

%% Extrude path to create surface

for ii=1:N
    [X(ii).manual, Y(ii).manual, Z(ii).manual] = extrude(base, path(ii).manual, 0);
    [X(ii).mag,    Y(ii).mag,    Z(ii).mag]    = extrude(base, path(ii).mag, 0);
    [X(ii).nomag,  Y(ii).nomag,  Z(ii).nomag]  = extrude(base, path(ii).nomag, 0);
end

%% Plot shaded colorbar for each trial
figure(1); clf(1); hold on; 

% plot extrusions and map forces to CData
for ii=1:N
    h(ii).manual = surf(X(ii).manual, Y(ii).manual, Z(ii).manual);
    h(ii).manual.CData = repmat(phantom(ii).manual.Fmag_smooth', [length(base), 1]);

    h(ii).mag = surf(X(ii).mag, Y(ii).mag, Z(ii).mag);
    h(ii).mag.CData = repmat(phantom(ii).mag.Fmag_smooth', [length(base), 1]);

    h(ii).nomag = surf(X(ii).nomag, Y(ii).nomag, Z(ii).nomag);
    h(ii).nomag.CData = repmat(phantom(ii).nomag.Fmag_smooth', [length(base), 1]);
end

% Add case labels
text_locations = [30*ones(size(y_locations)); y_locations; 1.05*width*ones(size(y_locations))]';
labels = {'Trial 1', 'Trial 2', 'Trial 3', 'Trial 4','Trial 1', 'Trial 2', 'Trial 3', 'Trial 4','Trial 1', 'Trial 2', 'Trial 3', 'Trial 4'};

% Find max force
temp = [phantom.manual];
forces = vertcat(temp.Fmag_smooth);
temp = [phantom.nomag];
forces = vertcat(forces, temp.Fmag_smooth);
temp = [phantom.mag];
forces = vertcat(forces, temp.Fmag_smooth);
max_force = max(forces);

% Colormap
load('forces_colormap.mat'); % loads forces_colormap (use colormapeditor to create/edit)
cmap = colormap(forces_colormap);

ax = gca;
ax.FontSizeMode = 'manual';
ax.FontSize = 14;
ax.XTick = 0:90:630;
ax.YTick = [];
xlabel('Angular Insertion Depth (deg)', 'FontSize',14);
xlim([10 630]);

set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
ylim([min(y_locations)-width/2, max(y_locations)+width/2]);

shading flat;

cbar = colorbar('Limits', [0 max_force]);
cbar.Label.String = 'Force (mN)';
cbar.Label.FontSize = 14';