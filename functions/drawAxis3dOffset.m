function drawAxis3dOffset(varargin)
% Modified from geom3d toolbox
% Trevor Bruns
% July 2019
%
%
%DRAWAXIS3D Draw a coordinate system and an origin
%
%   drawAxis3d
%	Adds 3 cylinders to the current axis, corresponding to the directions
%	of the 3 basis vectors Ox, Oy and Oz.
%	Ox vector is red, Oy vector is green, and Oz vector is blue.
%
%   drawAxis3d(origin)
%   Specifies the origin and used default axes.
%
%   drawAxis3d(origin, x_axis, y_axis, z_axis)
%   Specifies the origin and axes.
%
%   drawAxis3d(origin, x_axis, y_axis, z_axis, L, R)
%   Specifies the origin, axes, length L and the radius of the cylinders 
%   representing the different axes.
%
%   Example
%   drawAxis
%
%   figure;
%   drawAxis([1,1,1]);
%
%   See also
%   drawAxisCube
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2007-08-14,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2007 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.

% geometrical data
origin = [0 0 0];
v1 = [1 0 0];
v2 = [0 1 0];
v3 = [0 0 1];

% default parameters
L = 1;
r = L/10;

% extract parameters from input
if ~isempty(varargin)
	
    
    if length(varargin) <=3
        T = varargin{1};
        origin = T(1:3,4)';
        v1 = T(1:3,1)';
        v2 = T(1:3,2)';
        v3 = T(1:3,3)';
        if length(varargin) == 3
            L = varargin{2};
            r = varargin{3};
        end
    elseif length(varargin) > 3
        origin = varargin{1};
        v1 = varargin{2};
        v2 = varargin{3};
        v3 = varargin{4};

        if length(varargin) == 6
            L = varargin{5};
            r = varargin{6};
        end
    end
    
end



% draw 3 cylinders and a ball
hold on;
drawCylinder([origin origin+v1*L r], 16, 'facecolor', 'r', 'edgecolor', 'none');
drawCylinder([origin origin+v2*L r], 16, 'facecolor', 'g', 'edgecolor', 'none');
drawCylinder([origin origin+v3*L r], 16, 'facecolor', 'b', 'edgecolor', 'none');
drawSphere([origin 2*r], 'faceColor', 'black');

