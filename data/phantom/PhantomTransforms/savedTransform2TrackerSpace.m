function T_tracker_a = savedTransform2TrackerSpace(T_tracker_a)
    
    lps2ras = diag([-1,-1,1,1]);
    % Pull in the stored itk transform (stores it as an inverse)
    T_tracker_a = reshape(T_tracker_a,3,4);
    T_tracker_a = [T_tracker_a;[0,0,0,1]];
    
    % Similarity Transform to express transformation in lps space
    T_tracker_a = lps2ras*T_tracker_a*inv(lps2ras);
    T_tracker_a(1:3,4) = T_tracker_a(1:3,1:3)*T_tracker_a(1:3,4);
    
    % Flip translation vector
    T_tracker_a(1:3,4) = -T_tracker_a(1:3,4);
    
    % Doing all of this makes the saved transforms match what is output in
    % slicer and track!


end
