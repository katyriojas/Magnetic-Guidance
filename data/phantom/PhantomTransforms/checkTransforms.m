% Katy Riojas
% Created: 8/19/19
% Last Updated: 8/30/19
% Check Transformation Matrices

clear all; close all; clc;
ras2lps = diag([-1,-1,1,1]);

% Pull in desired transforms (move mag and ait to fixture 
% frame expressed in each of their frames)
% T_mag_fixture_goal = [1.00, 0.00, 0.00, 0.59;...
%                       0.00, 1.00, 0.00, 108.95;...
%                       0.00, 0.00, 1.00, -22.34;...
%                       0.00, 0.00, 0.00, 1.00];
                  
% T_ait_fixture_goal = [1.00, -0.02, 0.03, -0.00;...
%                       0.02, 1.00, 0.02, 6.08;...
%                       -0.04, -0.02, 1.00, -28.90;...
%                       0.00, 0.00, 0.00, 1.00];

fileID = fopen('T_ait_fixture.txt','r');
tform_ait = textscan(fileID, '%f, %f, %f, %f');
T_ait_fixture_goal = [tform_ait{1},tform_ait{2},tform_ait{3},tform_ait{4}];
fclose(fileID);

fileID = fopen('T_mag_fixture.txt','r');
tform_mag = textscan(fileID, '%f, %f, %f, %f');
T_mag_fixture_goal = [tform_mag{1},tform_mag{2},tform_mag{3},tform_mag{4}];
fclose(fileID);

% Pull in Transforms
basepath = strcat(pwd,'\ug-mea-saline-1.25\trial4-pre\');
errname = 'ug4-pre-';
toggleMagSave = 0;

load(strcat(basepath,'T_cochlea_tracker.mat'));
T_tracker_fixture = AffineTransform_double_3_3;

load(strcat(basepath,'T_ait_tracker.mat')); % Loads in affine transform_3_3
T_tracker_ait = AffineTransform_double_3_3;

load(strcat(basepath,'T_mag_tracker.mat'));
T_tracker_mag = AffineTransform_double_3_3;

% Convert Transform to Tracker Space
T_tracker_fixture = savedTransform2TrackerSpace(T_tracker_fixture);
T_tracker_ait = savedTransform2TrackerSpace(T_tracker_ait);
T_tracker_mag = savedTransform2TrackerSpace(T_tracker_mag);

T_tracker_ait_goal = T_tracker_fixture*inv(T_ait_fixture_goal);
T_tracker_mag_goal = T_tracker_fixture*inv(T_mag_fixture_goal); % v1 in python code

T_trackerait_trackeraitgoal = inv(T_tracker_ait)*T_tracker_ait_goal;
T_trackermag_trackermaggoal = inv(T_tracker_mag)*T_tracker_mag_goal;

Rait_tracked_goal = T_trackerait_trackeraitgoal(1:3,1:3);
Rmag_tracked_goal = T_trackermag_trackermaggoal(1:3,1:3);

ait_ang_err = vectorAngle3d([0,0,1],(Rait_tracked_goal*[0,0,1]')');
mag_ang_err = vectorAngle3d([0,0,1],(Rmag_tracked_goal*[0,0,1]')');

ait_ang_err_deg = rad2deg(ait_ang_err)
mag_ang_err_deg = rad2deg(mag_ang_err);

ait_tip_err = vecnorm(T_trackerait_trackeraitgoal(1:3,4));
mag_center_err = vecnorm(T_trackermag_trackermaggoal(1:3,4));

magcolor = 'b';
aitcolor = 'k';

figure(1);
subplot(1,2,1); grid on; hold on; title('Angular Offset');
xlabel('Trial'); ylabel('Angular Offset (axang rep) [deg]');
scatter(1,ait_ang_err_deg,aitcolor,'filled');
scatter(1,mag_ang_err_deg,magcolor,'filled');

subplot(1,2,2); grid on; hold on; title('Origin Offset');
xlabel('Trial'); ylabel('Origin Offset [mm]');
scatter(1,ait_tip_err,aitcolor,'filled');
scatter(1,mag_center_err,magcolor,'filled');
legend('NMAIT','Omnimag');

save(strcat(pwd,'\errors\ait_tip_err\',errname,'ait_tip_err.mat'),'ait_tip_err');
save(strcat(pwd,'\errors\ait_ang_err\',errname,'ait_ang_err_deg.mat'),'ait_ang_err_deg');

if toggleMagSave
    save(strcat(pwd,'\errors\mag_center_err\',errname,'mag_center_err.mat'),'mag_center_err');
    save(strcat(pwd,'\errors\mag_ang_err\',errname,'mag_ang_err_deg.mat'),'mag_ang_err_deg');
end

