function [pCochleaPPR, pEntryPPR, Fiducial_Markers, Extenders, voxel_size, ct_dims] = Read_PPR_File(PPR_FilePath)

if nargin == 0
    [PPR_filename, PPR_pathname] = uigetfile('.ppr','Select .ppr File');
    PPR_FilePath = fullfile(PPR_pathname, PPR_filename);
end

% Open file, and scan for lines containing markers
fileID = fopen(PPR_FilePath,'r');
M = [];
while ~feof(fileID)
    tline = fgets(fileID);
    if length(tline) >= 5
        if strcmp(tline(1:5), '[MARK')
           ScanningFids = 1;
           while ScanningFids == 1
               marker_row_txt = fgets(fileID);
               if length(marker_row_txt) > 10
                   M_tmp = textscan(marker_row_txt,'%f');
                   M_tmp = M_tmp{1};
                   M = [M, M_tmp]; %#ok<AGROW>
               else
                   ScanningFids = 0;
               end
           end           
        end
    end
end
fclose(fileID);

% Open file and scan for line trajectory data
fileID = fopen(PPR_FilePath,'r');
while ~feof(fileID)
    tline = fgets(fileID);
    if length(tline) >= 5
        if strcmp(tline(1:5), '[TRAJ')
           traj_row_txt = fgets(fileID);
           T4 = textscan(traj_row_txt,'%f');
        end
    end
end
fclose(fileID);

Traj = T4{1};
pCochleaPPR = Traj(1:3);
pEntryPPR = Traj(4:6);

%FiducialMarkers = zeros(3,3);
if nargout > 2
    if ~isempty(M) 
        Fiducial_Markers = M(1:3,:);
        Extenders = M(4:6,:);
    else
        warning('No marker data present in ppr file (not a problem if this is a preop scan).');
        Fiducial_Markers = zeros(3,3);
        Extenders = zeros(3,3);
    end
end

% Open file and scan for voxel size and ct dimensions
% NOTE: will only return info for first scan listed (e.g. preop/postop)
if nargout == 6
    fileID = fopen(PPR_FilePath,'r');
    while ~feof(fileID)
        tline = fgets(fileID);
        if length(tline) >= 15
            if strcmp(tline(6:15), 'IMAGE INFO')
               ct_info = cell2mat( textscan(tline(18:end),'%f') );
               ct_dims    = ct_info(1:3);
               voxel_size = ct_info(4:6);
               break;
            end
        end
    end
    fclose(fileID);
end