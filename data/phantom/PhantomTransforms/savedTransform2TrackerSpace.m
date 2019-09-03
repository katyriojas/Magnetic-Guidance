function T_a_tracker = savedTransform2TrackerSpace(T_a_tracker)
    
    lps2ras = diag([-1,-1,1,1]);
    % Pull in the stored itk transform (stores it as an inverse)
    T_a_tracker = reshape(T_a_tracker,3,4);
    T_a_tracker = [T_a_tracker;[0,0,0,1]];
    
    % Similarity Transform to express transformation in lps space
    T_a_tracker = lps2ras*T_a_tracker*inv(lps2ras);
    T_a_tracker(1:3,4) = T_a_tracker(1:3,1:3)*T_a_tracker(1:3,4);
    
    % Flip translation vector -- need to figure out why this is happening
    T_a_tracker(1:3,4) = -T_a_tracker(1:3,4);
    
    % Doing all of this makes the saved transforms match what is output in
    % slicer and track!


end
