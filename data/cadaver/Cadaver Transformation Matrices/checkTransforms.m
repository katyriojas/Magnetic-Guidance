% Katy Riojas
% Created: 8/19/19
% Last Updated: 8/30/19
% Check Transformation Matrices

% Note on 8/20 - it matters whether I calculated the error with the tracker
% positions or with the transformed relative positions -- look into this
% but how it is currently written it matches the python code

clear all; close all; clc;
ras2lps = diag([-1,-1,1,1]);

% T_mag_fixture_goal = [1.00, 0.00, 0.00, 2.68;...
%                       0.00, 1.00, 0.00, 112.27;...
%                       0.00, 0.00, 1.00, -29.10;...
%                       0.00, 0.00, 0.00, 1.00];
                  
% T_ait_fixture_goal = 0.90, -0.40, -0.16, 3.83
%                      0.28, 0.83, -0.47, 23.01
%                      0.32, 0.38, 0.87, -25.62
%                      0.00, 0.00, 0.00, 1.00];

fileID = fopen('T_ait_fixture.txt','r');
tform_ait = textscan(fileID, '%f, %f, %f, %f');
T_ait_fixture_goal = [tform_ait{1},tform_ait{2},tform_ait{3},tform_ait{4}];
fclose(fileID);

fileID = fopen('T_mag_fixture.txt','r');
tform_mag = textscan(fileID, '%f, %f, %f, %f');
T_mag_fixture_goal = [tform_mag{1},tform_mag{2},tform_mag{3},tform_mag{4}];
fclose(fileID);

% % Change inputs
basepath = strcat(pwd,'\g\trial3-pre\');
errname = 'g3-pre-';
toggleMagSave = 1;
%%
load(strcat(basepath,'T_cochlea_tracker.mat'));
T_fixture_tracker = AffineTransform_double_3_3;
load(strcat(basepath,'T_cochlea_tracker (2).mat'));
T_fixture_tracker2 = AffineTransform_double_3_3;
load(strcat(basepath,'T_cochlea_tracker (3).mat'));
T_fixture_tracker3 = AffineTransform_double_3_3;
T_fixture_tracker_avg = mean([T_fixture_tracker,T_fixture_tracker2,T_fixture_tracker3],2);

load(strcat(basepath,'T_ait_tracker.mat'));
T_ait_tracker = AffineTransform_double_3_3;
load(strcat(basepath,'T_ait_tracker (2).mat'));
T_ait_tracker2 = AffineTransform_double_3_3;
load(strcat(basepath,'T_ait_tracker (3).mat'));
T_ait_tracker3 = AffineTransform_double_3_3;
T_ait_tracker_avg = mean([T_ait_tracker,T_ait_tracker2,T_ait_tracker3],2);

load(strcat(basepath,'T_mag_tracker.mat'));
T_mag_tracker = AffineTransform_double_3_3;
load(strcat(basepath,'T_mag_tracker (2).mat'));
T_mag_tracker2 = AffineTransform_double_3_3;
load(strcat(basepath,'T_mag_tracker (3).mat'));
T_mag_tracker3 = AffineTransform_double_3_3;
T_mag_tracker_avg = mean([T_mag_tracker,T_mag_tracker2,T_mag_tracker3],2);

% Convert Transform to Tracker Space
T_fixture_tracker = savedTransform2TrackerSpace(T_fixture_tracker_avg);
T_ait_tracker = savedTransform2TrackerSpace(T_ait_tracker_avg);
T_mag_tracker = savedTransform2TrackerSpace(T_mag_tracker_avg);

T_ait_tracker_goal = T_fixture_tracker*inv(T_ait_fixture_goal);
T_mag_tracker_goal = T_fixture_tracker*inv(T_mag_fixture_goal); % v1 in python code

ait_tip_offset = T_ait_tracker(1:3,4) - T_ait_tracker_goal(1:3,4);
mag_center_offset = T_mag_tracker(1:3,4) - T_mag_tracker_goal(1:3,4);

ait_tip_offset_mag = vecnorm(ait_tip_offset);
mag_center_offset_mag = vecnorm(mag_center_offset);

R_ait_error = T_ait_tracker_goal(1:3,1:3)*T_ait_tracker(1:3,1:3)';
R_mag_error = T_mag_tracker_goal(1:3,1:3)*T_mag_tracker(1:3,1:3)';

axang_ait_err = rotm2axang(R_ait_error);
axang_mag_err = rotm2axang(R_mag_error);

ang_ait_err = rad2deg(axang_ait_err(4));
ang_mag_err = rad2deg(axang_mag_err(4));

magcolor = 'b';
aitcolor = 'k';

figure(1);
subplot(1,2,1); grid on; hold on; title('Angular Offset');
xlabel('Trial'); ylabel('Angular Offset (axang rep) [deg]');
scatter(1,ang_ait_err,aitcolor,'filled');
scatter(1,ang_mag_err,magcolor,'filled');

subplot(1,2,2); grid on; hold on; title('Origin Offset');
xlabel('Trial'); ylabel('Origin Offset [mm]');
scatter(1,ait_tip_offset_mag,aitcolor,'filled');
scatter(1,mag_center_offset_mag,magcolor,'filled');
legend('NMAIT','Omnimag');

save(strcat(pwd,'\errors\ait_tip_offset\',errname,'ait_mag.mat'),'ait_tip_offset_mag');
save(strcat(pwd,'\errors\R_ait_err\',errname,'ait_R.mat'),'R_ait_error');

if toggleMagSave
    save(strcat(pwd,'\errors\mag_cent_offset\',errname,'mag_mag.mat'),'mag_center_offset_mag');
    save(strcat(pwd,'\errors\R_mag_err\',errname,'mag_R.mat'),'R_mag_error');
end

