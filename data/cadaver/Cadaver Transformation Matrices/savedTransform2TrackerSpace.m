function T_tracker_a = savedTransform2TrackerSpace(T_tracker_a)
    
% Saving the transform node saves the inverse of the RAS coordinates and
% converts to LPS. 

    % Pull in the stored itk transform (stores it as an inverse)
    T_tracker_a = reshape(T_tracker_a,3,4);
    T_tracker_a = [T_tracker_a;[0,0,0,1]];
    
    % First we convert to RAS
    lps2ras = diag([-1,-1,1,1]);
    
    % Similarity Transform to express transformation in RAS slicer space
    T_tracker_a = lps2ras*T_tracker_a*inv(lps2ras);
    
    % Still need to handle the flip in order of rotation and translation
    T_tracker_a(1:3,4) = -T_tracker_a(1:3,1:3)*T_tracker_a(1:3,4);
    
    % Tinv = [RT -RT*p; 0,1];
    % Maybe it stores it as if it is translation and then rotation as
    % as opposed to rotation and then translation
    % Doing all of this makes the saved transforms match what is output in
    % slicer and track!
    
    
end
