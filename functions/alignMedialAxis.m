function [medial_axis_aligned, T] = alignMedialAxis(medial_axis, side, varargin)
%ALIGNMEDIALAXIS    Align a set of scala tympani medial axis points.
%   
%   This function takes in an ordered set of 3xN points that define the 
%   medial axis of the scala tympani (from base to apex) and aligns them to
%   the standard cochlear frame => 
%   Consensus panel on a cochlear coordinate system applicable in 
%   histological, physiological and radiological studies of the human
%   cochlea (Verbist et al., 2010)
%       
%   Right ear -> +Z axis towards helicotrema
%   Left  ear -> -Z axis towards helicotrema
%
%                       Y
%                       ^
%                      /|\
%                       |
%               ,aadPP""|""YYbawawawawawawawa     
%            ,aP"'      |       
%           aP'     ____|_____     
%          d"    ,adP"""|""""Yba,  
%        ,d'   ,dP"     |      `Yb,    
%       ,d'   ,d'    ,dP|"Yb,    `Y,     
%       d'    d'   ,d"  |   "b,   `Y,         
%       8     8    d'   +---------------------> X    
%       8     8    8     8    `8    8         
%       Y,    Y,   `b, ,aP     P    8      
%       `Y,   `Ya    """"     d'   ,P      
%        `8,    `Ya         ,8"   ,P'   
%         `Ya,    `Ya,,__,,d"'   ,P'      
%           `Ya,     `""""'     ,P'   
%             `"Ya,_          ,d"        
%                 ""YbaaaaaadP"           
%
%   MEDIAL_AXIS_ALIGNED = alignMedialAxis(MEDIAL_AXIS, SIDE, BASAL_PTS); 
%   returns the aligned points. BASAL_PTS are the indices of the medial
%   axis points to use when computing the central axis of the cochlea.
%
%   [MEDIAL_AXIS_ALIGNED, T] = alignMedialAxis(MEDIAL_AXIS, SIDE);
%   returns the aligned points the homogeneous transformation [4x4] which 
%   transforms points in the original frame to the aligned frame. 
%   (MEDIAL_AXIS_ALIGNED = T*MEDIAL_AXIS)
%
%
%   Trevor Bruns
%   June 2019

%% check inputs

% basal_pts = 5:40;
n_inputs = length(varargin);

if nargin < 2
    error('not enough input arguments')
end

if n_inputs == 1
    if length(varargin{1}) == 3
        insertion_axis = varargin{1};
    else
        basal_pts = varargin{2};
    end
elseif n_inputs == 2
    if length(varargin{1}) == 3
        insertion_axis = varargin{1};
        basal_pts = varargin{2};
    else
        basal_pts = varargin{1};
        insertion_axis = varargin{2};
    end
else
    error('too many input arguments')
end

if size(medial_axis,1) ~= 3
    error('medial_axis must be a [3xN] list of points')
end


%% Z-axis
% fit circle to basal turn
[center, z_axis, ~] = CircFit3D(medial_axis(:,basal_pts)');

% Make sure normal is oriented correctly (right ear -> towards helicotrema)
center_to_apex = normalizeVector3d((medial_axis(:,end) - center)')';
if xor( strcmpi(side, 'L'), (dot(z_axis, center_to_apex) < 0) )
    z_axis = -z_axis;
end

% create XY plane at the most basal point
xy_plane = createPlane((medial_axis(:,1))', z_axis');
% xy_plane = createPlane((medial_axis(:,1))', (medial_axis(:,1)+x_axis)', (medial_axis(:,1)+y_axis)');

% project apical end onto XY plane and fit circle to find center
apical_pts = projPointOnPlane(medial_axis(:,end-17:end)', xy_plane);
[origin, ~, ~] = CircFit3D(apical_pts);
    

%% X-axis
if ~exist('insertion_axis')
    % fit line to first N points
    num_straight_points = 15; % number of points to use for fitting line
    straight_points = medial_axis(:,1:num_straight_points);
    difference = bsxfun(@minus, straight_points', mean(straight_points, 2)');
    [~,~,V] = svd(difference, 0);
    insertion_axis = V(:,1); 
end

% project insertion axis onto XY plane
x_proj = projPointOnPlane((medial_axis(:,1)-insertion_axis)' , xy_plane)';
x_axis = x_proj - medial_axis(:,1);

%% Y-axis
y_axis = normalizeVector3d(cross(z_axis', x_axis'))';

%% Origin
% find intersection of Z-axis with the XY plane
% z_axis_line = createLine3d(center', (center + z_axis)');
% origin = intersectLinePlane(z_axis_line, xy_plane)';

%% Rotation matrix
% need R_aligned_orig, which rotates points in the original frame
% into the aligned frame. Our computed axes (i.e. basis vectors) are the 
% columns of R_orig_aligned. Thus, we simply take the transpose of this:
% R_aligned_original = inv(R_original_aligned) = R_original_aligned'
R_original_aligned = [x_axis, y_axis, z_axis];
R_original_aligned = quat2rotm(quatnormalize(rotm2quat(R_original_aligned))); % orthonormalize
R_aligned_original = R_original_aligned'; % now transpose

%% Align
% translate points by the origin offset and apply rotation
medial_axis_aligned = R_aligned_original * (medial_axis - repmat(origin,[1 length(medial_axis)]));

%% T
% A homogeneous transformation applies a rotation followed a translation.
% We have done the reverse. Thus we cannot simply create T by combining
% R_aligned_original and origin, e.g. T = [R, -origin; [0,0,0,1]]
%
% However, since we now have the aligned points, we can compute T via point
% registration:
if nargout == 2
    [R,t] = point_register(medial_axis,medial_axis_aligned);
    T = [R,t;[0,0,0,1]];
end

end