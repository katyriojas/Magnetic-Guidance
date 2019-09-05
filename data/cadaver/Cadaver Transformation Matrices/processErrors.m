% Katy Riojas
% 9/5/19
% Process transform error data
clear all; close all; clc;

ait_tip_dir = 'errors\ait_tip_offset\';
ait_tip_offset_files = dir(strcat(ait_tip_dir,'*.mat'));
ait_tip_offset_data  = cell(1, numel(ait_tip_offset_files));

for ii = 1:numel(ait_tip_offset_files)
  
    FileData1 = load(strcat(ait_tip_dir,ait_tip_offset_files(ii).name));
    ait_tip_offset_data{ii} = FileData1.ait_tip_offset_mag;
    
end

ait_tip_offset_mat = cat(1, ait_tip_offset_data{:}); 
maxAITTipErr = max(ait_tip_offset_mat);
fprintf('Cadaver Trials: Maximum ait 3d tip error %2.2f\n', maxAITTipErr);

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
fprintf('Cadaver Trials: Maximum magnet center error %2.2f\n', maxMagCenterErr);

%%

R_mag_dir = 'errors\R_mag_err\';
R_mag_files = dir(strcat(R_mag_dir,'*.mat'));
R_mag_data  = cell(1, numel(R_mag_files));

for ii = 1:numel(R_mag_files)
  
    FileData3 = load(strcat(R_mag_dir ,R_mag_files(ii).name));
    R_mag_data{ii} = FileData3.R_mag_error;
    mag_angles(ii,:) = rad2deg(rotm2eul(R_mag_data{ii}));
end

R_mag_data_mat= cat(1, mag_angles); %store along the rows
max_mag_Rerr = max(R_mag_data_mat(:));
fprintf('Cadaver Trials: Maximum euler angular offset of Magnet %2.2f deg\n',max_mag_Rerr);

%%
R_ait_dir = 'errors\R_ait_err\';
R_ait_files = dir(strcat(R_ait_dir,'*.mat'));
R_ait_data  = cell(1, numel(R_mag_files));

for ii = 1:numel(R_ait_files)
  
    FileData4 = load(strcat(R_ait_dir ,R_ait_files(ii).name));
    R_ait_data{ii} = FileData4.R_ait_error;
    ait_angles(ii,:) = rad2deg(rotm2eul(R_ait_data{ii}));
    
end

R_ait_data_mat= cat(1, ait_angles);
max_ait_Rerr = max(R_ait_data_mat(:));
fprintf('Cadaver Trials: Maximum euler angular offset of AIT %2.2f deg\n',max_ait_Rerr);