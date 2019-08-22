function markers = LoadFiducialMarkerLocations(tool_definition_filepath)
% LoadFiducialMarkerLocations - Retrieves fiducial markers locations from
% NDI tool definition file (.txt version)
%
%   markers = LoadFiducialMarkerLocations() asks the user to select a tool 
%   definition file and returns the marker positions as a [3xN] array
%   
%
% Example:
%   >> markers = LoadFiductialMarkerLocations();
%   >> markers = LoadFiductialMarkerLocations("C:\path\to\tool_definition.txt");
%
% Trevor Bruns
% June 2019

if nargin == 0
    [filename, pathname] = uigetfile('.txt','Select NDI Rigid Body File');
    tool_definition_filepath = fullfile(pathname, filename);
end

% open file
fileID = fopen(tool_definition_filepath,'r');

% scan until line before first marker position => "----------------"
keep_scanning = true;
while keep_scanning
    current_line = fgets(fileID); % read next line
    if strcmp(current_line(1), '-') % check if we have reached the line
        keep_scanning = false;
    end
    
    if feof(fileID) % reached end of file
        error('Reached end of file without match. Did you select the correct file?')
    end
end

% parse and store marker positions
markers = [];
while length(current_line) > 2 % continue until empty line after marker positions
    current_line = fgets(fileID); % read next line
    next_marker = cell2mat(textscan(current_line,'%f')); % convert to array
    markers = [markers, next_marker];
end

% close file
fclose(fileID);

end