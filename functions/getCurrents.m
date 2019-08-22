function [currents,condition_n] = getCurrents(tool_position,B_field)

    global M mu_0;
    tool_position_n = tool_position/norm(tool_position); %normalized tool position
    %Calculate Q matrix - maps from field applied to coil currents
    Q = (2*pi/mu_0)*(norm(tool_position))^3*(M^-1)...
                   *(3*tool_position_n...
                   *(tool_position_n') - 2*eye(3)); 
                 
    currents = Q*B_field; %assume putting 3 A through each coil
    
    if nargout > 1
%         [V,D] = eig(Q);
        condition_n = cond(Q);
    end
end