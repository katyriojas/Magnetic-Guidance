% Katy Riojas
% 9/10/19
% Process transform error data
clear all; close all; clc;

ait_tip_dir = 'errors\ait_tip_err\';
ait_tip_err_files = dir(strcat(ait_tip_dir,'*.mat'));
ait_tip_err_data  = cell(1, numel(ait_tip_err_files));

for ii = 1:numel(ait_tip_err_files)
  
    FileData1 = load(strcat(ait_tip_dir,ait_tip_err_files(ii).name));
    ait_tip_err_data{ii} = FileData1.ait_tip_err;
    
end

ait_tip_offset_mat = cat(1, ait_tip_err_data{:}); 
max_ait_tip_err = max(ait_tip_offset_mat);
fprintf('Cadaver Trials: Maximum AIT 3d tip error %2.2f\n', max_ait_tip_err);

%% Now Compute the magnet error
mag_cent_dir = 'errors\mag_center_err\';
mag_center_err_files = dir(strcat(mag_cent_dir,'*.mat'));
mag_cent_err_data = cell(1, numel(mag_center_err_files));

for ii = 1:numel(mag_center_err_files)
  
    FileData2 = load(strcat(mag_cent_dir ,mag_center_err_files(ii).name));
    mag_cent_err_data{ii} = FileData2.mag_center_err;
    
end

mag_cent_offset_mat = cat(1, mag_cent_err_data{:});
max_mag_center_err = max(mag_cent_offset_mat);
fprintf('Cadaver Trials: Maximum magnet center error %2.2f\n', max_mag_center_err);

%% Angular Insertion Tool Error
ait_ang_err_dir = 'errors\ait_ang_err\';
ait_ang_err_files = dir(strcat(ait_ang_err_dir,'*.mat'));
ait_ang_err_data = cell(1, numel(ait_ang_err_files));

for ii = 1:numel(ait_ang_err_files)
  
    FileData3 = load(strcat(ait_ang_err_dir ,ait_ang_err_files(ii).name));
    ait_ang_err_data{ii} = FileData3.ait_ang_err_deg;
    
end

ait_ang_err_mat= cat(1, ait_ang_err_data{:});
max_ait_ang_err = max(ait_ang_err_mat(:));
fprintf('Cadaver Trials: Maximum angular offset of AIT %2.2f deg\n',max_ait_ang_err);

%% Angular Magnet Error
mag_ang_err_dir = 'errors\mag_ang_err\';
mag_ang_err_files = dir(strcat(mag_ang_err_dir,'*.mat'));
mag_ang_err_data = cell(1, numel(mag_ang_err_files));

for ii = 1:numel(mag_ang_err_files)
  
    FileData4 = load(strcat(mag_ang_err_dir ,mag_ang_err_files(ii).name));
    mag_ang_err_data{ii} = FileData4.mag_ang_err_deg;
    
end

mag_ang_err_mat= cat(1, mag_ang_err_data{:});
max_mag_ang_err = max(mag_ang_err_mat(:));
fprintf('Cadaver Trials: Maximum angular offset of Omnimagnet is %2.2f deg\n',max_mag_ang_err);