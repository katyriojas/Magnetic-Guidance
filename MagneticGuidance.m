%% Creates Magnetic Guidance Plan
% 
% Trevor Bruns & Katy Riojas
% June/July 2019
%
% ----------
% --INPUTS--
% ----------
% ST medial axis  => .txt file of points along medial axis (from IMPROVISE)
% insertion trajectory => .ppr file from CIP with target/entry/marker pts
% cochlea fixture tool definition => .txt NDI tool definition file
% 
% -----------
% --OUTPUTS--
% -----------
% T_ait_fixture.txt => target pose of AIT in cochlea fixture frame
% T_mag_fixture.txt => target pose of Omnimagnet in cochlea fixture frame
% magnetic_guidance_plan.yaml => coil currents and insertion depths
%
% Notes:
% - requires geom3d, geom2d, robotics, and aerospace toolboxes
%   - geom3d = matlab.addons.toolbox.installToolbox('geom3d.mltbx')
%   - geom2d = matlab.addons.toolbox.installToolbox('geom2d.mltbx')

clear all; clc; close all;
load('path.mat');

addpath('stls','functions','omnimag parameters','rig_ppr_medial_path',...
    'preoperative plans');

%% User-Specified Parameters

% Toggles
align_st_with_mag = false;
useLW = true;
flip_magnet_polarity = true; %true equal south out the tip
ramp_field = true;

% Left or Right Ear
% side = 'R';
side = 'L';

if side == 'L'
    path(3,:) = -path(3,:);
end

% Specify medial_axis interpolation step size and start/end depths
interp_step = 0.04; % [mm]
start_depth = 8; % [mm]
end_depth = 27; %[mm]

insertion_speed = 1.25; %[mm/s]
% [mm] distance between center of the ST and Omnimagnet's surface
st_offset_from_mag_surface = 16; 

% position [mm] of scala tympani (center) in the Omnimagnet's coordinate frame
% 84 is hardcoded in here for half of the width of the Omnimagnet
t_mag_st = [0; 84 + st_offset_from_mag_surface; 0];
        
% Desired magnitude of magnetic field (this assumes constant along trajectory)
Bmag = 78*10^-3; % [T] (experimentally: 104mT yields ~75mT)
scale1 = 5/(Bmag*1000); % assuming plus or minus 5 mT
Bmag_start = 70*10^-3; %[mT]
Bmag_end = 85*10^-3; % [mT]
ramp_start = 16; % [mm] (along insertion depth)
ramp_end = 22.5; %[mm] (along insertion depth)
% Select which points to use for fitting the basal plane
% basal_pts = 15:40; % medial axis points to use when fitting basal plane
basal_pts = 40:80;

% STL paths
ait_path = 'AIT_7-17-19.STL';
mag_path = 'omnimag_7-17-19.STL';
fixture_path = 'cochlea_fixturev4.STL';

current_scaling_x = 1.31;
current_scaling_y = 1.74*1.06;
current_scaling_z = 1.38*1.082;

%% Load STLs

ait_stl = stlRead(ait_path);
mag_stl = stlRead(mag_path);
fixture_stl = stlRead(fixture_path);


%% Load insertion trajectory (.ppr), scala tympani medial axis (.txt), and cochlea fixture tool definition file (.txt)

% Specify file locations

% PPR file location
% ppr_filepath = 'rig_ppr_medial_path\phantom1_preopPlan.ppr';
ppr_filepath = 'rig_ppr_medial_path\tbone4.ppr';

% Medial axis file location
% medial_axis_filepath = 'rig_ppr_medial_path\phantom1_medial_axis_ct.txt';
medial_axis_filepath = 'rig_ppr_medial_path\MedialAxis_Tbone4.txt';

% Fixture file location
fixture_filepath = strcat('C:\Users\riojaske\Documents\magsteer\',...
'PreopPlanning\rig_ppr_medial_path\CochleaFixtureRef_2019-7-17.txt');

% Load data from ppr file (uses RPI coordinate frame)
[target_rpi, entry_rpi, markers_rpi, ~, voxel_size, ct_dims] =...
    Read_PPR_File(ppr_filepath);

% Load medial axis points
medial_axis_lpi = LoadMedialAxis(side, medial_axis_filepath);
              
% Load fiducial marker locations from cochlea fixture tool definition
markers_fixture = LoadFiducialMarkerLocations(fixture_filepath);


%% Transform CT data to LPS coordinate frame

% PPR data: convert from RPI to LPS
current_frame = 'RPI';
desired_frame = 'LPS';

% First convert markers to LPS
num_markers = size(markers_rpi,2);
markers_lps = zeros(3,num_markers);

for ii = 1:num_markers
    
    markers_lps(:,ii) = ...
        transform_between_coordinate_systems_Matlab(markers_rpi(1:3,ii),...
        ct_dims, voxel_size, current_frame, desired_frame);
    
end

% Now convert the target and entry points to LPS
target_lps = transform_between_coordinate_systems_Matlab(target_rpi,...
    ct_dims, voxel_size, current_frame, desired_frame);
entry_lps = transform_between_coordinate_systems_Matlab(entry_rpi,...
    ct_dims, voxel_size, current_frame, desired_frame);
target_vector_lps = normalizeVector3d(target_lps' - entry_lps')';

% Convert Medial Axis Points from LPI (Image Voxels) to LPS
medial_axis_lpi_mm = bsxfun(@times, medial_axis_lpi, voxel_size); % voxels to mm
medial_axis_lps_mm = zeros(size(medial_axis_lpi)); % initialize LPS medial axis

% Compute medial axis points in LPS in mm
for ii = 1:length(medial_axis_lpi)
    
    medial_axis_lps_mm(1:3,ii) =...
        transform_between_coordinate_systems_Matlab(medial_axis_lpi_mm(1:3,ii),...
        ct_dims, voxel_size, 'LPI', 'LPS');
    
end


%% Compute registration from ct frame to cochlea fixture frame => T_fixture_lps

% Markers from fixture frame come from NDI rig body txt file and
% markers_lps comes from CT space
[T_fixture_lps, FRE] = ...
    point_register_without_correspondence(markers_lps, markers_fixture);

if FRE > 1.0
    warning('FRE is large (%.2f mm). Did you select the correct files?',FRE)
end


%% Align medial axis to standard scala tympani coordinate frame
[medial_axis_st, T_st_lps] = alignMedialAxis(medial_axis_lps_mm,...
    side,basal_pts, target_vector_lps);

% Also compute T_fixture_st (needed later)
T_fixture_st = T_fixture_lps * inv(T_st_lps);
T_st_fixture = inv(T_fixture_st);

%% Interpolate and smooth medial axis to add points and ensure they are evenly spaced
smooth_span = 0.1;
medial_axis_smoothed = interp_and_smooth(medial_axis_st,interp_step,smooth_span);
lateral_wall_smoothed = interp_and_smooth(path,interp_step,smooth_span);

if ~useLW
    path_diff = sqrt(sum(diff(medial_axis_smoothed')'.^2, 1));
else
    path_diff = sqrt(sum(diff(lateral_wall_smoothed')'.^2,1));
end

path_cumulative = [0, cumsum(path_diff)]; % cumulative length along path
path_length_mm = path_cumulative(end); % total path length

if ramp_field
    
    ramp_start_index = find(path_cumulative>=ramp_start,1);
    ramp_end_index = find(path_cumulative>=ramp_end,1);
    Bmag_ramp = linspace(Bmag_start,Bmag_start,length(path_cumulative));
    Bmag_ramp(ramp_start_index:ramp_end_index) = ...
        linspace(Bmag_start,Bmag_end,length(ramp_start_index:ramp_end_index));
    Bmag_ramp(ramp_end_index:end) = Bmag_end;
    
end

%% Compute magnetic field vectors in cochlea frame

% compute alpha => angle of normal vector w.r.t. x-axis
% Note: ignoring z-component since vector should lie mostly in XY plane
medial_axis_gradient = gradient(medial_axis_smoothed);
lateral_wall_gradient = gradient(lateral_wall_smoothed);
alpha = bsxfun(@atan2d, medial_axis_gradient(2,:),...
            medial_axis_gradient(1,:)) - 90;
alpha2 = bsxfun(@atan2d, lateral_wall_gradient(2,:),...
            lateral_wall_gradient(1,:)) - 90;

% normal vectors at each point (negative since they point inward)
if useLW~=1
    mag_vectors_st = [-cosd(alpha); -sind(alpha); zeros(size(alpha))];
else
    mag_vectors_st = [-cosd(alpha2); -sind(alpha2); zeros(size(alpha2))];
end

if flip_magnet_polarity
   mag_vectors_st = -mag_vectors_st; 
end

%% Transform magnetic field vectors into Omnimagnet coordinate frame
% Need both the vector origins (i.e. medial axis points) and vectors in 
% Omnimagnet frame

% The cochlea fixture and the Omnimagnet are fixed parallel to the table 
% surface, so we can only rotate about zhat_fixture in order to best 
% align yhat_st with yhat_mag

 a = projPointOnPlane((T_fixture_st(1:3,1:3)*[0;1;0])',...
     createPlane([0 0 0], [0 0 1]))';

% assume y axis of of fixture is initially aligned with mag before rotation
b = [0;1;0]; 

if align_st_with_mag
    
    R_theta_z = createRotationVector3d(a', [0,1,0]); % compute angular offset
    
else
    
    R_theta_z = eye(4);
    
end

% Rotation about Z-axis that aligns ST Y-axis with the fixture Y-axis

t_theta_z = R_theta_z*T_fixture_st;
t_theta_z = t_theta_z(1:3,4);

T_mag_fixture = R_theta_z;
T_mag_fixture(1:3,4) =  t_mag_st - t_theta_z;

% T_mag_st => transformation from scala tympani frame to Omnimagnet frame
T_mag_st = T_mag_fixture * T_fixture_st;

% Transform the magnetic field vectors into Omnimagnet frame
mag_vectors_mag = T_mag_st(1:3,1:3) * mag_vectors_st; % can just rotate the field vectors
medial_axis_mag = transformPoint3d(medial_axis_smoothed', T_mag_st)';

% Compute for Lateral Wall array
lateral_wall_mag = transformPoint3d(lateral_wall_smoothed',T_mag_st)';

%% Compute omnimagnet coil currents

initialize_omnimag_vars;

if useLW 
    endcurrents = length(lateral_wall_mag);
    mag_pnts = lateral_wall_mag;
    currents = zeros(size(lateral_wall_mag));
else
    endcurrents = length(medial_axis_mag);
    mag_pnts = medial_axis_mag;
    currents = zeros(size(medial_axis_mag));
end

for ii = 1:endcurrents
    
    if ~ramp_field
        currents(:,ii) = getCurrents(mag_pnts(:,ii)/1000,...
        mag_vectors_mag(:,ii)*Bmag);
    else
        currents(:,ii) = getCurrents(mag_pnts(:,ii)/1000,...
        mag_vectors_mag(:,ii)*Bmag_ramp(ii));
    end
    
end

currents_scaled = zeros(size(currents));
% Scale the currents if necessary
currents_scaled(1,:) = current_scaling_x*currents(1,:);
currents_scaled(2,:) = current_scaling_y*currents(2,:);
currents_scaled(3,:) = current_scaling_z*currents(3,:);

%% Trim to desired start/end depths

mag_start = find(path_cumulative > start_depth, 1);
mag_end   = find(path_cumulative > end_depth,   1);

planned_pnts = mag_start:mag_end;
spacing = (end_depth - start_depth)/(length(planned_pnts)-1);
insertion_depth = start_depth:spacing:end_depth;


%% Find transform for AIT target pose in fixture coordinate frame
% T_fixture_ait = T_fixture_lps * (T_st_lps)^-1 * T_st_ait
%                    ?have          ?have           ?need

% First determine R_st_ait => rotation of the AIT target pose in the 
% scala tympani frame to do this we will find the basis vectors (axes)
% of the AIT frame in ST coordinates Z-axis (axis of insertion)
target_vector_st = T_st_lps(1:3,1:3)*target_vector_lps;

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

% Now we just need the translation component before assembling T_st_ait
% Assume AIT frame origin is the first trajectory/medial axis point
t_st_ait = medial_axis_smoothed(:,1); 
T_st_ait = [R_st_ait, t_st_ait; [0,0,0,1]];

% Compute T_fixture_ait
T_fixture_ait = T_fixture_lps * inv(T_st_lps) * T_st_ait;
T_ait_fixture = inv(T_fixture_ait);

%% Compute transform for cochlea fixture frame in Omnimagnet frame

T_fixture_mag = T_fixture_lps * inv(T_st_lps) * inv(T_mag_st);


%% Transform STLs into Omnimagnet frame

% AIT
T_mag_ait = inv(T_fixture_mag) * T_fixture_ait;
ait_stl_magframe = ait_stl;
ait_stl_magframe.vertices = transformPoint3d(ait_stl.vertices, T_mag_ait);

% Cochlea fixture
fixture_stl_magframe = fixture_stl;
fixture_stl_magframe.vertices = transformPoint3d(fixture_stl_magframe.vertices, inv(T_fixture_mag));


%% Estimate coil temperatures
start_temps = [25;25;25];
time_step = interp_step/insertion_speed; % [s] time between each point
coil_temps = omnimagnetCoilTempEstimate(currents_scaled(:,planned_pnts), time_step, start_temps);

% plot
% Omnimagnet coil temperatures
figure(4); clf(4)
grid on; hold on; axis equal; xlabel('x'); ylabel('y');
plot(path_cumulative(mag_start:mag_end), coil_temps(1,:),...
     path_cumulative(mag_start:mag_end), coil_temps(2,:),...
     path_cumulative(mag_start:mag_end), coil_temps(3,:));
ylim([min(coil_temps(:))-2, max(coil_temps(:))+2]);

for ii = 1:3
    [max_temp, max_ind] = max(coil_temps(ii,:));
    text(path_cumulative(max_ind + mag_start), max_temp, strcat(' \leftarrow ', sprintf(['%.1f' char(176) 'C'], max_temp)))
end

legend('Coil_x','Coil_y','Coil_z'); 
xlabel('Distance Along Cochlea Path (mm)');
ylabel('Temperature (C)');
title('Omnimagnet Coil Temperatures')

%% Plot

% Omnimagnet coil currents
figure(1); clf(1); 
hax=axes;
grid on; hold on;
plot(path_cumulative,currents_scaled(1,:),...
     path_cumulative,currents_scaled(2,:),...
     path_cumulative,currents_scaled(3,:));

% Add vertical lines for start/end depths
line([path_cumulative(mag_start) path_cumulative(mag_start)],...
    get(hax,'YLim'),'Color','green','LineStyle','--','LineWidth',1.5)
line([path_cumulative(mag_end) path_cumulative(mag_end)],...
    get(hax,'YLim'),'Color','red','LineStyle','--','LineWidth',1.5)

legend('Ix','Iy','Iz','start point','end point'); 
xlim([0 ceil(path_length_mm)]);
xlabel('Distance Along Cochlea Path (mm)');
ylabel('Current (A)');
[max_z, max_z_ind] = max(currents_scaled(3,:));
[min_z, min_z_ind] = min(currents_scaled(3,:));
text(path_cumulative(max_z_ind),max_z, strcat(' \leftarrow ',...
    sprintf('%.1f A', max_z)))
text(path_cumulative(min_z_ind),min_z, strcat(' \leftarrow ',...
    sprintf('%.1f A', min_z)))
title('Omnimagnet Coil Currents')

% 3D, scala tympani frame
figure(2); clf(2)
hold on; axis equal; xlabel('x'); ylabel('y'); zlabel('z');
drawPolyline3d(medial_axis_st','--r');
drawPolyline3d(medial_axis_smoothed', 'b');
drawPolyline3d(path','--g');
drawPolyline3d(lateral_wall_smoothed', 'k');
mag_on = mag_start:mag_end;
if ~useLW
    mag_off = 1:length(medial_axis_smoothed);
else
    mag_off = 1:length(medial_axis_smoothed);
end

mag_off(mag_on) = [];

if ~useLW
    drawFieldVectors(medial_axis_smoothed,mag_vectors_st,mag_off,mag_on)
else
    drawFieldVectors(lateral_wall_smoothed,mag_vectors_st,mag_off,mag_on)
end
    
drawAxis3d(0.5, 0.025);
drawAxis3dOffset(T_st_ait, 0.5, 0.025);
drawPlane3d([0,0,0, 1,0,0, 0,1,0], 'FaceColor', 'm', 'FaceAlpha', 0.1)
drawVector3d(t_st_ait', target_vector_st', 'LineWidth', 2.5, 'Color', 'm');
view([-100 35]);

% % 2D, scala tympani frame
% figure(5); clf(5)
% hold on; axis equal; xlabel('x'); ylabel('y'); zlabel('z');
% drawPolyline(medial_axis_st(1:2,:)', 'k');
% drawPolyline(medial_axis_smoothed(1:2,:)', '--b');
% drawVector(medial_axis_smoothed(1:2,1:2:end)', mag_vectors_st(1:2,1:2:end)','c');

% 3D, Omnimagnet frame
figure(3); clf(3); grid on;
hold on; axis equal; xlabel('x'); ylabel('y'); zlabel('z');

drawPoint3d(lateral_wall_mag','b');
drawPolyline3d(medial_axis_mag','k');

if ~useLW
    drawVector3d(medial_axis_mag(:,1:5:end)', mag_vectors_mag(:,1:5:end)','k');
else
    drawVector3d(lateral_wall_mag(:,1:5:end)', mag_vectors_mag(:,1:5:end)','k');
end
drawAxis3dOffset(T_mag_st, 10, 0.4);
drawCuboid([0,0,0, 178, 168, 172], 'FaceColor', [1 0.5 0], 'FaceAlpha', 0.15);
patch('Faces',mag_stl.faces, 'Vertices', mag_stl.vertices,...
    'FaceColor', [.7 0.7 .7], 'EdgeColor', 'none');
drawAxis3d(90,0.75);
patch('Faces',ait_stl_magframe.faces, 'Vertices',...
    ait_stl_magframe.vertices, 'FaceColor', [0 1 0.5], 'EdgeColor', 'none');
drawAxis3dOffset(T_mag_ait, 10, 0.4);
patch('Faces',fixture_stl_magframe.faces, 'Vertices',...
    fixture_stl_magframe.vertices, 'FaceColor', [0 0.5 0.7], 'EdgeColor',...
    'none','FaceAlpha', 0.3);
drawAxis3dOffset(T_mag_fixture, 12, 0.4);

% Aesthetics
view([123 11]);
set(gcf, 'renderer', 'opengl')
p_light = [150, 1200, 350];
light('Position', p_light, 'Style', 'local')
material dull

%% Export plan to .yaml file for ROS
newFolderName = strcat('preoperative plans/PreopPlan_',...
    datestr(now,'yyyy-mm-dd_HH-MM'));
mkdir(newFolderName);

% Save figures
savefig(figure(1),strcat(newFolderName,'\Fig1'));
savefig(figure(2),strcat(newFolderName,'\Fig2'));
savefig(figure(3),strcat(newFolderName,'\Fig3'));
savefig(figure(4),strcat(newFolderName,'\Fig4'));

yamlFileID = fopen(strcat(newFolderName,'/trajectory', '.yaml'),'w');
export_data2ros(1:length(planned_pnts),...
                insertion_depth,...
                currents_scaled(:,planned_pnts),yamlFileID);
            
%% Export T_fixture_ait.txt => target pose of AIT in cochlea fixture frame
fileID = fopen(strcat(newFolderName,'/T_ait_fixture','.txt'), 'wt');
fprintf(fileID, '%2.2f, %2.2f, %2.2f, %2.2f\n', T_ait_fixture(1,:));
fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T_ait_fixture(2,:));
fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T_ait_fixture(3,:));
fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T_ait_fixture(4,:));
fclose(fileID);


%% Export T_fixture_mag.txt => target pose of Omnimagnet in cochlea fixture frame
fileID = fopen(strcat(newFolderName,'/T_mag_fixture', '.txt'), 'wt');
fprintf(fileID, '%2.2f, %2.2f, %2.2f, %2.2f\n', T_mag_fixture(1,:));
fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T_mag_fixture(2,:));
fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T_mag_fixture(3,:));
fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T_mag_fixture(4,:));
fclose(fileID);

if export_data
    % Export plan to .yaml file for ROS

    newFolderName = strcat('preoperative plans/PreopPlan_',...
        datestr(now,'yyyy-mm-dd_HH-MM'));
    mkdir(newFolderName);

    % Save figures
    savefig(figure(1),strcat(newFolderName,'\Fig1'));
    savefig(figure(2),strcat(newFolderName,'\Fig2'));
    savefig(figure(3),strcat(newFolderName,'\Fig3'));
    savefig(figure(4),strcat(newFolderName,'\Fig4'));

    yamlFileID = fopen(strcat(newFolderName,'/trajectory', '.yaml'),'w');
    export_data2ros(1:length(planned_pnts),...
                    insertion_depth,...
                    currents_scaled(:,planned_pnts),yamlFileID);


    % Export T_fixture_ait.txt => target pose of AIT in cochlea fixture frame
    fileID = fopen(strcat(newFolderName,'/T_ait_fixture','.txt'), 'wt');
    fprintf(fileID, '%2.2f, %2.2f, %2.2f, %2.2f\n', T_ait_fixture(1,:));
    fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T_ait_fixture(2,:));
    fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T_ait_fixture(3,:));
    fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T_ait_fixture(4,:));
    fclose(fileID);

    % Export T_fixture_mag.txt => target pose of Omnimagnet in cochlea fixture frame
    fileID = fopen(strcat(newFolderName,'/T_mag_fixture', '.txt'), 'wt');
    fprintf(fileID, '%2.2f, %2.2f, %2.2f, %2.2f\n', T_mag_fixture(1,:));
    fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T_mag_fixture(2,:));
    fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T_mag_fixture(3,:));
    fprintf(fileID,  '%2.2f, %2.2f, %2.2f, %2.2f\n', T_mag_fixture(4,:));
    fclose(fileID);

    % Export Matlab workspace
    save(strcat(newFolderName,'/workspace.mat'));
end

% Commented out - used the code below to generate preop plan figure
% figure(5); grid on; hold on;
% % ylim([66 90]);
% xlabel('Insertion Depth (mm)');
% ylabel('Bmag (mT)');
% plot(insertion_depth,Bmag_ramp(planned_pnts)*1000,'k','LineWidth',1);
% figure();
% grid on; hold on;
% plot(insertion_depth,currents_scaled(1,planned_pnts),...
%      insertion_depth,currents_scaled(2,planned_pnts),...
%      insertion_depth,currents_scaled(3,planned_pnts));
% xlabel('Insertion Depth (mm)'); ylabel('Currents (A)');
% legend('Ix','Iy','Iz');