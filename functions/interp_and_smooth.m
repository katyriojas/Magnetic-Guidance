function medial_axis_smoothed = interp_and_smooth(medial_axis_st,interp_step,smooth_span)

medial_axis_diff = sqrt(sum(diff(medial_axis_st')'.^2, 1));
medial_axis_cumsum = [0, cumsum(medial_axis_diff)]; % cumulative length along path
interp_pts = 0:interp_step:medial_axis_cumsum(end);
medial_axis_interp = [interp1(medial_axis_cumsum, medial_axis_st(1,:), interp_pts);
                      interp1(medial_axis_cumsum, medial_axis_st(2,:), interp_pts);
                      interp1(medial_axis_cumsum, medial_axis_st(3,:), interp_pts)];

medial_axis_smoothed = [ smooth(medial_axis_interp(1,:), smooth_span, 'loess')'; ...
                         smooth(medial_axis_interp(2,:), smooth_span, 'loess')'; ...
                         smooth(medial_axis_interp(3,:), smooth_span, 'loess')' ];
                         
end