function offset_poly = offsetPolyline(poly, dist)
% offsets a polyline (2xN) by a specified distance


% compute alpha => angle of normal vector w.r.t. x-axis
poly_gradient = gradient(poly);
alpha = bsxfun(@atan2d, poly_gradient(2,:), poly_gradient(1,:)) - 90;


% normal vectors at each point 
normals = [cosd(alpha); sind(alpha)];

% apply offset
offset_poly = poly + dist*normals;

end