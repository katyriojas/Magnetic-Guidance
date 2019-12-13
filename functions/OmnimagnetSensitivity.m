function [pos_error_mag, pos_error_ang, ang_error_mag, ang_error_ang] = OmnimagnetSensitivity(p_nom, b_nom)
% Sensitivity analysis of dipole-field source pose
% Jake J. Abbott

mu0 = 4*pi*1e-7;

p_hat = normalizeVector3d(p_nom')';

% compute nominal m, using nominal B at nominal p
m_nom = (2*pi/mu0) * norm(p_nom)^3 * (3*(p_hat*p_hat') - 2*eye(3)) * b_nom;


%%
% First, let's look at the effect of position error of the dipole source
% Create 1-mm radius sphere of p vectors surrounding nominal p.
% Evaluate field at each p location, and then calculate vector error,
% magnitude error (as a percentage), and angular error (in degrees) of
% nominal assumption.

n = 30;
rho = 0.001;
theta = linspace(0,2*pi,n);
phi = linspace(0,pi,n);
p = zeros(3,n*n);
b = zeros(3,n*n);
vector_error = zeros(3,n*n);
magnitude_error = zeros(1,n*n);
angle_error = zeros(1,n*n);
index = 0;
for i = 1:n
    for j = 1:n
        index = index + 1;
        p(:,index) = p_nom + [rho*sin(phi(j))*cos(theta(i)); rho*sin(phi(j))*sin(theta(i)); rho*cos(phi(j))];
        b(:,index) = mu0/(4*pi*norm(p(:,index))^3)*(3/(norm(p(:,index))^2)*p(:,index)*p(:,index)' - eye(3))*m_nom;
        vector_error(:,index) = b(:,index) - b_nom;
        magnitude_error(:,index) = abs((norm(b(:,index)) - norm(b_nom)))/norm(b(:,index))*100;
        angle_error(:,index) = acos(b(:,index)'*b_nom/(norm(b_nom)*norm(b(:,index))))*180/pi;
    end
end

% Report conservative error bounds for position error in placement of
% dipole source
max_magnitude_error_percentage_from_position_error = max(magnitude_error);
max_angle_error_degrees_from_position_error        = max(angle_error);


%%
% Next, let's look at the effect of angular errors of the dipole source
% Create 1-degree error in m vector (in two dimensions) from nominal m
% First, make orthogonal vectors from m

seed_vect = rand(3,1);
vect1 = seed_vect - (seed_vect'*m_nom/norm(m_nom))*m_nom/norm(m_nom);
vect1 = vect1/norm(vect1);
vect2 = cross(m_nom,vect1);
n = 30;
psi = 1*pi/180; % This is the 1-degree magnitude
theta = linspace(0,pi/2,n);
phi = linspace(-psi,psi,3);
m = zeros(3,n*3);
b = zeros(3,n*3);
vector_error = zeros(3,n*3);
magnitude_error = zeros(1,n*3);
angle_error = zeros(1,n*3);
index = 0;
for i = 1:n
    for j = 1:3
        index = index + 1;
        k = cos(theta(i))*vect1 + sin(theta(i))*vect2;
        skew = [0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];
        R = expm(phi(j)*skew);
        m(:,index) = R*m_nom;
        b(:,index) = mu0/(4*pi*norm(p_nom)^3)*(3/(norm(p_nom)^2)*(p_nom*p_nom') - eye(3))*m(:,index);
        vector_error(:,index) = b(:,index) - b_nom;
        magnitude_error(:,index) = abs((norm(b(:,index)) - norm(b_nom)))/norm(b(:,index))*100;
        angle_error(:,index) = acos(b(:,index)'*b_nom/(norm(b_nom)*norm(b(:,index))))*180/pi;
    end
end

% Report conservative error bounds for position error in placement of
% dipole source
max_magnitude_error_percentage_from_angular_error = max(magnitude_error);
max_angle_error_degrees_from_angular_error	      = max(angle_error);



%% Assemble outputs

pos_error_mag = max_magnitude_error_percentage_from_position_error;
pos_error_ang = max_angle_error_degrees_from_position_error;
ang_error_mag = max_magnitude_error_percentage_from_angular_error;
ang_error_ang = max_angle_error_degrees_from_angular_error;

end