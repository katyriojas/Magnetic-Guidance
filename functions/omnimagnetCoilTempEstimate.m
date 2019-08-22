function coil_temps = omnimagnetCoilTempEstimate(coil_currents, time_step, start_temps)

if nargin < 3
    start_temps = [25;25;25]; % [C]
end

% Omnimagnet properties
R = [3.5, 3.8, 4.0]; % [ohms] resistances of inner (x), middle (y), and outer (z) coils 
area = 1.60967e-6; % [m^2] 16ga square copper wire from MWS
resistivity = 16.78e-9; % resistivity of copper
density = 8960; % [kg/m^3] density of copper
C = 385.0; % [J/(kg*C)] specific heat of copper
l = R.*(area/resistivity); % [m] lengths of coil windings
volume = l.*area; % [m^3] volume of coil windings
mass = volume.*density; % [kg] mass of coil windings

time = linspace(0, time_step*length(coil_currents), length(coil_currents));
% time = 0:time_step:length(coil_currents);

cu_alpha = 0.004041; % copper temperature coefficient
coil_temps = zeros(size(coil_currents));
coil_temps(:,1) = start_temps;

for ii_coil = 1:3
    Q_cumulative = 0; % [J] cumulative energy converted to heat

    for ii_time=1:(length(time)-1)
        % compute energy generated during next time step
        P_curr = coil_currents(ii_coil,ii_time)^2 * ...
            (R(ii_coil)*(1+cu_alpha*(coil_temps(ii_coil,ii_time)-...
            start_temps(ii_coil))));

        % 'integrate'
        Q_cumulative = Q_cumulative + P_curr*time_step;

        % compute current temperature
        dt = Q_cumulative / (mass(ii_coil)*C);
        coil_temps(ii_coil,ii_time+1) = start_temps(ii_coil) + dt;
        
    end
end

end