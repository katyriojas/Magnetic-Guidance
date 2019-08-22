function lh_nomag = addData2Plot(data_nomag,alpha,colors)
    
    % unguided, no smoothing
    lh_nomag(1) = plot(data_nomag.depth_insertion, data_nomag.Fx, 'Color', [colors(1,:),0.3*alpha]);
    lh_nomag(2) = plot(data_nomag.depth_insertion, data_nomag.Fy, 'Color', [colors(2,:),0.3*alpha]);
    lh_nomag(3) = plot(data_nomag.depth_insertion, data_nomag.Fz, 'Color', [colors(3,:),0.3*alpha]);

    % unguided, smoothed
    lh_nomag(4) = plot(data_nomag.depth_insertion, data_nomag.Fx_smooth, 'Color', [colors(1,:),0.99*alpha]);
    lh_nomag(5) = plot(data_nomag.depth_insertion, data_nomag.Fy_smooth, 'Color', [colors(2,:),0.99*alpha]);
    lh_nomag(6) = plot(data_nomag.depth_insertion, data_nomag.Fz_smooth, 'Color', [colors(3,:),0.99*alpha]);
    
    set(lh_nomag(1:3), 'LineWidth', 0.8)
    set(lh_nomag(4:6), 'LineWidth', 1)

end