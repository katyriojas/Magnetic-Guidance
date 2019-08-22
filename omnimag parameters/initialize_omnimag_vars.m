%initialize Magnet things
global M mu_0 %we will need these in getCurrents and getField

% Need the linear transformation to map coil currents to dipole
load('OptimalSqCore.mat'); %this will load in the matrix of values.
coreMoment = optimalGeometry.CoreMoment;
solenoidMoments = optimalGeometry.SolenoidMoment;
summedMoments = coreMoment + solenoidMoments; %[Am}
Lout = 25.4*optimalGeometry.OuterLength; %[m]
Wout = 25.4*optimalGeometry.OuterWidth;  %[m]
Tout = 25.4*optimalGeometry.OuterThickness;  %[m]
W2Tout = Wout + 2*Tout;
W2T = W2Tout;

Lmid = 25.4*optimalGeometry.MiddleThickness;  %[m]
Wmid = 25.4*optimalGeometry.MiddleWidth;  %[m]
Tmid = 25.4*optimalGeometry.MiddleThickness;  %[m]
W2Tmid = Wmid + 2*Tmid;

% Assuming 16 AWG square self-bonding copper wire: 
% https://mwswire.com/wp-content/uploads/2017/08/Microsquare-Magnet
% -Wire-MWS-Wire-Industries.pdf
wire_cross_section = 2495; %[mil^2]
wire_cross_section = wire_cross_section*6.4516e-10; %[m^2]

%M is a linear transformation that maps the three applied coil currents
%to the dipole moment m. [(A*m^2)/A]
M = summedMoments/wire_cross_section;
M = diag(M);

%Magnetic Constant for Vacuum Permeabiltiy 
%[m*kg*s^-2*A^-2] [also known as 4pi_e-7 [T*m*A^-1]
mu_0 = 4*pi*10^-7;