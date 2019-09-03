% Katy Riojas
% Created: 8/19/19
% Last Updated: 8/30/19
% Check Transformation Matrices

% Note on 8/20 - it matters whether I calculated the error with the tracker
% positions or with the transformed relative positions -- look into this
% but how it is currently written it matches the python code

clear all; close all; clc;
ras2lps = diag([-1,-1,1,1]);

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
basepath = strcat(pwd,'\ug-ea-saline-1.25\trial1-post\');
errname = 'ug1-ea-post-';
toggleMagSave = 0;

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

%%

ait_tip_dir = 'errors\ait_tip_offset\';
ait_tip_offset_files = dir(strcat(ait_tip_dir,'*.mat'));
ait_tip_offset_data  = cell(1, numel(ait_tip_offset_files));

for ii = 1:numel(ait_tip_offset_files)
  
    FileData1 = load(strcat(ait_tip_dir,ait_tip_offset_files(ii).name));
    ait_tip_offset_data{ii} = FileData1.ait_tip_offset_mag;
    
end

ait_tip_offset_mat = cat(1, ait_tip_offset_data{:}); 
maxAITTipErr = max(ait_tip_offset_mat);

%%

mag_cent_dir = 'errors\mag_cent_offset\';
mag_cent_offset_files = dir(strcat(mag_cent_dir,'*.mat'));
mag_cent_offset_data  = cell(1, numel(mag_cent_offset_files ));

for ii = 1:numel(mag_cent_offset_files)
  
    FileData2 = load(strcat(mag_cent_dir ,mag_cent_offset_files(ii).name));
    mag_cent_offset_data{ii} = FileData2.mag_center_offset_mag;
    
end

mag_cent_offset_mat = cat(1, mag_cent_offset_data{:});
maxMagCenterErr = max(mag_cent_offset_mat);

%%

R_mag_dir = 'errors\R_mag_err\';
R_mag_files = dir(strcat(R_mag_dir,'*.mat'));
R_mag_data  = cell(1, numel(R_mag_files));

for ii = 1:numel(R_mag_files)
  
    FileData3 = load(strcat(R_mag_dir ,R_mag_files(ii).name));
    R_mag_data{ii} = FileData3.R_mag_error;
    
end

R_mag_data_mat= cat(3, R_mag_data{:});


%%
R_ait_dir = 'errors\R_ait_err\';
R_ait_files = dir(strcat(R_ait_dir,'*.mat'));
R_ait_data  = cell(1, numel(R_mag_files));

for ii = 1:numel(R_ait_files)
  
    FileData4 = load(strcat(R_ait_dir ,R_ait_files(ii).name));
    R_ait_data{ii} = FileData4.R_ait_error;
    
end

R_ait_data_mat= cat(3, R_ait_data{:});