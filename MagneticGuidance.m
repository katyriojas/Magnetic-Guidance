%% Creates Magnetic Guidance Plan
% 
% Trevor Bruns & Katy Riojas
% June-September 2019
%
% 
% -----------
% -- NOTES --
% -----------
% - requires geom3d, geom2d, robotics, and aerospace toolboxes
%   - geom3d = matlab.addons.toolbox.installToolbox('geom3d.mltbx')
%   - geom2d = matlab.addons.toolbox.installToolbox('geom2d.mltbx')
% - requires 'Border-less tight subplot' from File Exchange
%   - https://www.mathworks.com/matlabcentral/fileexchange/68326-border-less-tight-subplot-auto-refresh
% - verify all toggles/flags are correct
% - check file paths
% - check that basal_pts make sense for the medial_axis file used
%
%
% ------------
% -- INPUTS --
% ------------
% ST medial axis  => .txt file of points along medial axis (from IMPROVISE)
% insertion trajectory => .ppr file from CIP with target/entry/marker pts
% cochlea fixture tool definition => .txt NDI tool definition file
% 
% -------------
% -- OUTPUTS --
% -------------
% T_ait_fixture.txt => target pose of AIT in cochlea fixture frame
% T.mag_fixture.txt => target pose of Omnimagnet in cochlea fixture frame
% magnetic_guidance_plan.yaml => coil currents and insertion depths
%

% clear all; clc; close all;
addpath('stls','functions','omnimag parameters','rig_ppr_medial_path', 'preoperative plans');


%% User-Specified Parameters

% Toggles
params.align_st_with_mag = false;   % true => align y-axes of ST & magnet; false => keep fixture parallel to magnet
params.useLW = true;                % true => use lateral wall path for generating field vectors
params.flip_magnet_polarity = true; % true => south out the tip
params.ramp_field = true;           % true => use ramp function for magnetic field magnitudes
params.export_data = false;         % true => create new folder and export generated data

params.ppr_filename = 'tbone4.ppr'; % PPR file location
% params.ppr_filename = 'phantom1_preopPlan.ppr'; % PPR file location


params.medial_axis_filename = 'MedialAxis_Tbone4.txt'; % Medial axis file location
% params.medial_axis_filename = 'phantom1_medial_axis_ct.txt'; % Medial axis file location

params.basal_pts = 40:80; % medial axis points to use when fitting basal plane


params.side = 'L'; % Left or Right Ear

load('path.mat'); % loads 'path' -> lateral wall path
if strcmp(params.side,'L')
    path(3,:) = -path(3,:);
end

% Specify medial_axis interpolation step size and start/end depths
params.interp_step = 0.04; % [mm]
params.start_depth = 8; % [mm]
params.end_depth = 27; %[mm]


params.insertion_speed = 1.25; % [mm/s] linear insertion speed

params.st_offset_from_mag_surface = 16; % [mm] distance between ST center and Omnimagnet's surface


% Position [mm] of scala tympani (center) in the Omnimagnet's coordinate frame
% 84 is hardcoded in here for half of the width of the Omnimagnet
params.t_mag_st = [0; 84 + params.st_offset_from_mag_surface; 0];
        

% Desired magnetic field strength parameters    
params.Bmag_const = 78e-3; % [T]  (used if params.ramp_field = false)

params.Bmag_start = 70e-3; % [T] (used if params.ramp_field = true)
params.ramp_start = 16;    % [mm] (along insertion depth)
params.Bmag_end   = 85e-3; % [T]
params.ramp_end   = 22.5;  % [mm] (along insertion depth)


% STL paths
params.ait_stl_filename     = 'AIT_7-17-19.STL';
params.mag_stl_filename     = 'omnimag_7-17-19.STL';
params.fixture_stl_filename = 'cochlea_fixturev4.STL';

% Fixture rigid body file location
params.fixture_filename = 'CochleaFixtureRef_2019-7-17.txt';

% Experimentally measured amounts to scale coil currents (NOTE: location dependent!!)
params.current_scaling.x = 1.31;
params.current_scaling.y = 1.84;
params.current_scaling.z = 1.49;


%% Load STLs

stl.ait = stlRead(params.ait_stl_filename);
stl.omnimag = stlRead(params.mag_stl_filename);
stl.fixture = stlRead(params.fixture_stl_filename);


%% Load insertion trajectory (.ppr), scala tympani medial axis (.txt), and cochlea fixture tool definition file (.txt)

% Load data from ppr file (uses RPI coordinate frame)
[ppr.target_rpi, ppr.entry_rpi, ppr.markers_rpi, ~, ppr.voxel_size, ppr.ct_dims] = Read_PPR_File(params.ppr_filename);

% Load medial axis points
medial_axis.lpi = LoadMedialAxis(params.side, params.medial_axis_filename);
              
% Load fiducial marker locations from cochlea fixture tool definition
markers.fixture = LoadFiducialMarkerLocations(params.fixture_filename);


%% Transform CT data to LPS coordinate frame

% PPR data: convert from RPI to LPS
current_frame = 'RPI';
desired_frame = 'LPS';

% First convert markers to LPS
num_markers = size(ppr.markers_rpi,2);
markers.lps = zeros(3,num_markers);

for ii = 1:num_markers
    
    markers.lps(:,ii) = ...
        transform_between_coordinate_systems_Matlab(ppr.markers_rpi(1:3,ii), ppr.ct_dims, ppr.voxel_size, current_frame, desired_frame);
    
end

% Now convert the target and entry points to LPS
target_lps = transform_between_coordinate_systems_Matlab(ppr.target_rpi, ppr.ct_dims, ppr.voxel_size, current_frame, desired_frame);
entry_lps  = transform_between_coordinate_systems_Matlab(ppr.entry_rpi,  ppr.ct_dims, ppr.voxel_size, current_frame, desired_frame);
target_vector_lps = normalizeVector3d(target_lps' - entry_lps')';

% Convert Medial Axis Points from LPI (Image Voxels) to LPS
medial_axis.lpi_mm = bsxfun(@times, medial_axis.lpi, ppr.voxel_size); % voxels to mm
medial_axis.lps_mm = zeros(size(medial_axis.lpi)); % initialize LPS medial axis

% Compute medial axis points in LPS in mm
for ii = 1:length(medial_axis.lpi)
    
    medial_axis.lps_mm(1:3,ii) = ...
        transform_between_coordinate_systems_Matlab(medial_axis.lpi_mm(1:3,ii), ppr.ct_dims, ppr.voxel_size, 'LPI', 'LPS');
    
end


%% Compute registration from ct frame to cochlea fixture frame => T.fixture_lps

% Markers from fixture frame come from NDI rig body txt file and
% markers.lps comes from CT space
[T.fixture_lps, FRE] =  point_register_without_correspondence(markers.lps, markers.fixture);

if FRE > 1.0
    warning('FRE is large (%.2f mm). Did you select the correct files?',FRE)
end


%% Align medial axis to standard scala tympani coordinate frame
[medial_axis.st, T.st_lps] = alignMedialAxis(medial_axis.lps_mm, params.side, params.basal_pts, target_vector_lps);

% Also compute T.fixture_st (needed later)
T.fixture_st = T.fixture_lps * inv(T.st_lps);
% T.st_fixture = inv(T.fixture_st);


%% Interpolate and smooth medial axis to add points and ensure they are evenly spaced
params.smooth_span = 0.1;
medial_axis.st_smoothed = interp_and_smooth(medial_axis.st, params.interp_step, params.smooth_span);
lateral_wall.smoothed = interp_and_smooth(path, params.interp_step, params.smooth_span);

if params.useLW
    path_diff = sqrt(sum(diff(lateral_wall.smoothed')'.^2, 1));
else
    path_diff = sqrt(sum(diff(medial_axis.st_smoothed')'.^2, 1));
end

insertion_depth = [0, cumsum(path_diff)]; % cumulative length along path
path_length_mm = insertion_depth(end); % total path length

if params.ramp_field
    ramp_start_index = find(insertion_depth>=params.ramp_start,1);
    ramp_end_index = find(insertion_depth>=params.ramp_end,1);
    Bmag_ramp = params.Bmag_start * ones(size(insertion_depth));
    Bmag_ramp(ramp_start_index:ramp_end_index) = ...
        linspace(params.Bmag_start, params.Bmag_end, length(ramp_start_index:ramp_end_index));
    Bmag_ramp(ramp_end_index:end) = params.Bmag_end;
end


%% Compute magnetic field vectors in cochlea frame

if params.useLW
    % compute alpha => angle of normal vector w.r.t. x-axis
    % Note: ignoring z-component since vector should lie mostly in XY plane
    lateral_wall.gradient = gradient(lateral_wall.smoothed);
    alpha = bsxfun( @atan2d, lateral_wall.gradient(2,:), lateral_wall.gradient(1,:) ) - 90;

    % normal vectors at each point (negative since they point inward)
    mag.vectors_st = [-cosd(alpha); -sind(alpha); zeros(size(alpha))];
else
    % compute alpha => angle of normal vector w.r.t. x-axis
    % Note: ignoring z-component since vector should lie mostly in XY plane
    medial_axis.st_smoothed_gradient = gradient(medial_axis.st_smoothed);
    alpha = bsxfun( @atan2d, medial_axis.st_smoothed_gradient(2,:), medial_axis.st_smoothed_gradient(1,:) ) - 90;
    
    % normal vectors at each point (negative since they point inward)
    mag.vectors_st = [-cosd(alpha); -sind(alpha); zeros(size(alpha))];
end

if params.flip_magnet_polarity
   mag.vectors_st = -mag.vectors_st; 
end


%% Transform magnetic field vectors into Omnimagnet coordinate frame
% Need both the vector origins (i.e. medial axis points) and vectors in 
% Omnimagnet frame

% The cochlea fixture and the Omnimagnet are fixed parallel to the table 
% surface, so we can only rotate about zhat_fixture in order to best 
% align yhat_st with yhat_mag

a = projPointOnPlane((T.fixture_st(1:3,1:3)*[0;1;0])', createPlane([0 0 0], [0 0 1]))';

% assume y axis of of fixture is initially aligned with mag before rotation
b = [0;1;0]; 

if params.align_st_with_mag
    R_theta_z = createRotationVector3d(a', [0,1,0]); % compute angular offset
else
    R_theta_z = eye(4);
end

% Rotation about Z-axis that aligns ST Y-axis with the fixture Y-axis

t_theta_z = R_theta_z*T.fixture_st;
t_theta_z = t_theta_z(1:3,4);

T.mag_fixture = R_theta_z;
T.mag_fixture(1:3,4) =  params.t_mag_st - t_theta_z;

% T.mag_st => transformation from scala tympani frame to Omnimagnet frame
T.mag_st = T.mag_fixture * T.fixture_st;

% Transform the magnetic field vectors into Omnimagnet frame
mag.vectors_mag = T.mag_st(1:3,1:3) * mag.vectors_st; % can just rotate the field vectors
medial_axis.mag = transformPoint3d(medial_axis.st_smoothed', T.mag_st)';
lateral_wall.mag = transformPoint3d(lateral_wall.smoothed',T.mag_st)';

if params.useLW
    mag.pnts = lateral_wall.mag;
else
    mag.pnts = medial_axis.mag;
end

%% Compute Omnimagnet coil currents

% Find indices where magnet is on/off
mag.i_start = find(insertion_depth > params.start_depth, 1);
mag.i_end   = find(insertion_depth > params.end_depth,   1);
mag.on = mag.i_start:mag.i_end;
mag.off = 1:length(insertion_depth);
mag.off(mag.on) = [];

% Desired field magnitude as function of insertion depth
Bmag = params.ramp_field*Bmag_ramp + (~params.ramp_field)*ones(size(insertion_depth))*params.Bmag_const;
Bmag(mag.off) = 0;

% compute coil currents
initialize_omnimag_vars;
currents.ideal = zeros(size(mag.vectors_mag));
for ii = 1:length(insertion_depth)
    currents.ideal(:,ii) = getCurrents(mag.pnts(:,ii)/1000, mag.vectors_mag(:,ii)*Bmag(ii));
end

% Scale currents using experimentally determined scale factors
currents.scaled = zeros(size(currents.ideal));
currents.scaled(1,:) = params.current_scaling.x * currents.ideal(1,:);
currents.scaled(2,:) = params.current_scaling.y * currents.ideal(2,:);
currents.scaled(3,:) = params.current_scaling.z * currents.ideal(3,:);


%% Find transform for AIT target pose in fixture coordinate frame
% T.fixture_ait = T.fixture_lps * (T.st_lps)^-1 * T.st_ait
%                    ^have          ^have           ^need

% First determine R_st_ait => rotation of the AIT target pose in the 
% scala tympani frame to do this we will find the basis vectors (axes)
% of the AIT frame in ST coordinates Z-axis (axis of insertion)
target_vector_st = T.st_lps(1:3,1:3)*target_vector_lps;

% Z axis of the AIT target pose in ST frame; AIT inserts along -Z axis
z_axis_st_ait = -target_vector_st; 

% Y-axis (electrode array curl axis)
xy_plane_ait = createPlane([0,0,0], z_axis_st_ait'); % create plane normal to z_axis_st_ait
y_proj = projPointOnPlane([0,1,0], xy_plane_ait); % curl points along Y axis of AIT frame
y_axis_st_ait = normalizeVector3d(y_proj)';

% X-axis
x_axis_st_ait = cross(y_axis_st_ait, z_axis_st_ait); % orthogonal to y & z axes

% R_st_ait
R_st_ait = [x_axis_st_ait, y_axis_st_ait, z_axis_st_ait];

% Now we just need the translation component before assembling T.st_ait
% Assume AIT frame origin is the first trajectory/medial axis point
t_st_ait = medial_axis.st_smoothed(:,1); 
T.st_ait = [R_st_ait, t_st_ait; [0,0,0,1]];

% Compute T.fixture_ait
T.fixture_ait = T.fixture_lps * inv(T.st_lps) * T.st_ait;
T.ait_fixture = inv(T.fixture_ait);


%% Compute transform for cochlea fixture frame in Omnimagnet frame

T.fixture_mag = T.fixture_lps * inv(T.st_lps) * inv(T.mag_st);


%% Transform STLs into Omnimagnet frame

% AIT
T.mag_ait = inv(T.fixture_mag) * T.fixture_ait;
stl.ait_magframe = stl.ait;
stl.ait_magframe.vertices = transformPoint3d(stl.ait.vertices, T.mag_ait);

% Cochlea fixture
stl.fixture_magframe = stl.fixture;
stl.fixture_magframe.vertices = transformPoint3d(stl.fixture_magframe.vertices, inv(T.fixture_mag));


%% Estimate coil temperatures
start_temps = [25;25;25];
time_step = params.interp_step/params.insertion_speed; % [s] time between each point
coil_temps = omnimagnetCoilTempEstimate(currents.scaled(:,mag.on), time_step, start_temps);


%% Plot

if exist('h_fig','var')
    if isvalid(h_fig)
        close(h_fig)
    end
end

h_fig = figure;
h_fig.WindowState = "maximized";

%%%
% Desired Bmag
h_ax(1) = subplot_er(3,2,1);
grid on;
plot(insertion_depth, Bmag*1e3, 'Color','k', 'LineWidth',1.5)
ylim([0,100])
xlabel('Distance Along Cochlea Path (mm)');
ylabel('||B|| [mT]');
title('Desired Magnetic Field', 'FontSize',13)


%%%
% Omnimagnet coil currents
h_ax(2) = subplot_er(3,2,3);
grid on; hold on;

plot(insertion_depth, currents.scaled(1,:), 'LineWidth',1.5)
plot(insertion_depth, currents.scaled(2,:), 'LineWidth',1.5)
plot(insertion_depth, currents.scaled(3,:), 'LineWidth',1.5)

% Add vertical lines for start/end depths
line([insertion_depth(mag.i_start) insertion_depth(mag.i_start)], get(h_ax(2),'YLim'),'Color','green', 'LineStyle','--', 'LineWidth',1.5)
line([insertion_depth(mag.i_end)   insertion_depth(mag.i_end)],   get(h_ax(2),'YLim'),'Color','red',   'LineStyle','--', 'LineWidth',1.5)

legend('Ix','Iy','Iz', 'Location','sw'); 
xlim([0 ceil(path_length_mm)]);
xlabel('Distance Along Cochlea Path (mm)');
ylabel('Current (A)');
[I.max_z, I.max_z_ind] = max(currents.scaled(3,mag.on));
[I.min_z, I.min_z_ind] = min(currents.scaled(3,mag.on));
[I.max_y, I.max_y_ind] = max(currents.scaled(2,mag.on));
[I.min_y, I.min_y_ind] = min(currents.scaled(2,mag.on));
[I.max_x, I.max_x_ind] = max(currents.scaled(1,mag.on));
[I.min_X, I.min_X_ind] = min(currents.scaled(1,mag.on));
text(insertion_depth(I.max_z_ind+mag.i_start),I.max_z, strcat('\leftarrow ', sprintf('%.1f A', I.max_z)))
text(insertion_depth(I.min_z_ind+mag.i_start),I.min_z, strcat('\leftarrow ', sprintf('%.1f A', I.min_z)))
text(insertion_depth(I.max_y_ind+mag.i_start),I.max_y, strcat('\leftarrow ', sprintf('%.1f A', I.max_y)))
text(insertion_depth(I.min_y_ind+mag.i_start),I.min_y, strcat('\leftarrow ', sprintf('%.1f A', I.min_y)))
text(insertion_depth(I.max_x_ind+mag.i_start),I.max_x, strcat('\leftarrow ', sprintf('%.1f A', I.max_x)))
text(insertion_depth(I.min_X_ind+mag.i_start),I.min_X, strcat('\leftarrow ', sprintf('%.1f A', I.min_X)))
title('Omnimagnet Coil Currents', 'FontSize',13)

%%%
% Omnimagnet coil temperatures
h_ax(3) = subplot_er(3,2,5);
grid on; hold on;
plot(insertion_depth(mag.on), coil_temps(1,:), 'LineWidth',1.5)
plot(insertion_depth(mag.on), coil_temps(2,:), 'LineWidth',1.5)
plot(insertion_depth(mag.on), coil_temps(3,:), 'LineWidth',1.5);
ylim([min(coil_temps(:))-2, max(coil_temps(:))+2]);

for ii = 1:3
    [max_temp, max_ind] = max(coil_temps(ii,:));
    text(insertion_depth(max_ind + mag.i_start), max_temp, strcat(' \leftarrow ', sprintf(['%.1f' char(176) 'C'], max_temp)))
end

% Add vertical lines for start/end depths
line([insertion_depth(mag.i_start) insertion_depth(mag.i_start)], get(h_ax(3),'YLim'),'Color','green', 'LineStyle','--', 'LineWidth',1.5)
line([insertion_depth(mag.i_end)   insertion_depth(mag.i_end)],   get(h_ax(3),'YLim'),'Color','red',   'LineStyle','--', 'LineWidth',1.5)

legend('T_x','T_y','T_z', 'Location','nw'); 
xlabel('Distance Along Cochlea Path (mm)');
ylabel(['Temperature (' char(176) 'C)']);
title('Omnimagnet Coil Temperatures', 'FontSize',13)

linkaxes(h_ax, 'x');

%%%
% 3D, scala tympani frame
subplot_er(3,2,2);
hold on; axis equal; xlabel('x'); ylabel('y'); zlabel('z');
drawPolyline3d(medial_axis.st','--r');
drawPolyline3d(medial_axis.st_smoothed', 'b');
drawPolyline3d(path','--g');
drawPolyline3d(lateral_wall.smoothed', 'k');


if params.useLW
    drawFieldVectors(lateral_wall.smoothed,mag.vectors_st,mag.off,mag.on)
else
    drawFieldVectors(medial_axis.st_smoothed,mag.vectors_st,mag.off,mag.on)
end
    
drawAxis3d(0.5, 0.025);
drawAxis3dOffset(T.st_ait, 0.5, 0.025);
drawPlane3d([0,0,0, 1,0,0, 0,1,0], 'FaceColor', 'm', 'FaceAlpha', 0.1)
drawVector3d(t_st_ait', target_vector_st', 'LineWidth', 2.5, 'Color', 'm');
view(2)


% % 2D, scala tympani frame
% figure(5); clf(5)
% hold on; axis equal; xlabel('x'); ylabel('y'); zlabel('z');
% drawPolyline(medial_axis.st(1:2,:)', 'k');
% drawPolyline(medial_axis.st_smoothed(1:2,:)', '--b');
% drawVector(medial_axis.st_smoothed(1:2,1:2:end)', mag.vectors_st(1:2,1:2:end)','c');


% 3D, Omnimagnet frame
subplot_er(3,2,[4,6]);
hold on; axis equal; xlabel('x'); ylabel('y'); zlabel('z');

drawPoint3d(lateral_wall.mag','b');
drawPolyline3d(medial_axis.mag','k');

if ~params.useLW
    drawVector3d(medial_axis.mag(:,1:5:end)',  mag.vectors_mag(:,1:5:end)','k');
else
    drawVector3d(lateral_wall.mag(:,1:5:end)', mag.vectors_mag(:,1:5:end)','k');
end
drawAxis3dOffset(T.mag_st, 10, 0.4);
drawCuboid([0,0,0, 178, 168, 172], 'FaceColor',[1 0.5 0], 'FaceAlpha',0.15);
patch('Faces',stl.omnimag.faces, 'Vertices',stl.omnimag.vertices, 'FaceColor',[.7 0.7 .7], 'EdgeColor','none');
drawAxis3d(90,0.75);
patch('Faces',stl.ait_magframe.faces, 'Vertices',stl.ait_magframe.vertices, 'FaceColor',[0 1 0.5], 'EdgeColor','none');
drawAxis3dOffset(T.mag_ait, 10, 0.4);
patch('Faces',stl.fixture_magframe.faces, 'Vertices',stl.fixture_magframe.vertices, 'FaceColor',[0 0.5 0.7], 'EdgeColor','none', 'FaceAlpha',0.3);
drawAxis3dOffset(T.mag_fixture, 12, 0.4);

% Aesthetics
view([123 11]);
set(gcf, 'renderer', 'opengl')
p_light = [150, 1200, 350];
light('Position', p_light, 'Style', 'local')
material dull


%% Export data

if params.export_data

    % Create new folder
    newFolderName = strcat('preoperative plans/PreopPlan_', datestr(now,'yyyy-mm-dd_HH-MM'));
    mkdir(newFolderName);

    % Save figures
    savefig(h_fig,strcat(newFolderName,'\plan'));
    
    % Export Matlab workspace
    delete(h_fig) % temporarily delete figure before exporting workspace since it is already saved
    save(strcat(newFolderName,'/workspace.mat'));
    openfig(strcat(newFolderName, '\plan.fig')); % re-open figure

    % Export T_fixture_ait.txt => target pose of AIT in cochlea fixture frame
    fileID = fopen(strcat(newFolderName,'/T_ait_fixture','.txt'), 'wt');
    fprintf(fileID, '%2.2f, %2.2f, %2.2f, %2.2f\n', T.ait_fixture(1,:));
    fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T.ait_fixture(2,:));
    fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T.ait_fixture(3,:));
    fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T.ait_fixture(4,:));
    fclose(fileID);


    % Export T_fixture_mag.txt => target pose of Omnimagnet in cochlea fixture frame
    fileID = fopen(strcat(newFolderName,'/T_mag_fixture', '.txt'), 'wt');
    fprintf(fileID, '%2.2f, %2.2f, %2.2f, %2.2f\n', T.mag_fixture(1,:));
    fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T.mag_fixture(2,:));
    fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T.mag_fixture(3,:));
    fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T.mag_fixture(4,:));
    fclose(fileID);

    % Export plan to .yaml file for ROS
    yamlFileID = fopen(strcat(newFolderName,'/trajectory', '.yaml'),'w');
    export_data2ros(insertion_depth(mag.on), currents.scaled(:,mag.on), yamlFileID);
    
end

%% Generate plan plot
figure(22); clf(22); 
xlimit = 27.05;
subplot(2,1,1); grid on; hold on;
set(gca,'FontSize',9,'FontName','Times');

plot(insertion_depth, currents.scaled(1,:), 'LineWidth',1.5,'Color','r')
plot(insertion_depth, currents.scaled(2,:), 'LineWidth',1.5,'Color','g')
plot(insertion_depth, currents.scaled(3,:), 'LineWidth',1.5,'Color','b')

legend('Ix','Iy','Iz', 'Location','nw','FontSize',9,'FontName','Times'); 
xlim([0,xlimit]);
set(gca,'xticklabel',{[]});
yticks([-50,-25,0,25,50]); %[A]
ylabel('Current (A)','FontSize',9,'FontName','Times','FontWeight','bold');

subplot(2,1,2); grid on; hold on;
set(gca,'FontSize',9,'FontName','Times');
plot(insertion_depth, Bmag*1e3, 'Color','k','LineWidth',1.5);
scatter(params.start_depth,params.Bmag_start*1000,5,'filled','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(params.ramp_start,params.Bmag_start*1000,5,'filled','MarkerFaceColor','k','MarkerEdgeColor','k');
scatter(params.ramp_end,params.Bmag_end*1000,5,'filled','MarkerFaceColor','k','MarkerEdgeColor','k');
ylim([0,100]);
yticks([0,25,50,75,100]); %[mT]
xlim([0,xlimit]); % [mm]
xlabel('Insertion Depth (mm)','FontSize',8,'FontName','Times','FontWeight','bold');
ylabel('||b|| (mT)','FontSize',8,'FontName','Times','FontWeight','bold');

h_Bmag_figure = gcf; 
h_Bmag_figure.PaperUnits = 'inches';
h_Bmag_figure.PaperSize = [8.5 11];
h_Bmag_figure.PaperPosition = [0 0 3.8 1.7];
saveas(h_Bmag_figure,'saved figures\currents_Bmag_plan.pdf');