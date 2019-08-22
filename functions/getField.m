function b_field = getField(tool_position,currents)

    global M mu_0;
    %normalized tool position
    tool_position_n = tool_position/norm(tool_position); 
    %Calculate Q matrix - maps from field applied to coil currents
    Q = (2*pi/mu_0)*(norm(tool_position))^3*(M^-1)...
                   *(3*tool_position_n...
                   *(tool_position_n') - 2*eye(3)); 

    b_field= Q^-1*currents; %assume putting 3 A through each coil

end