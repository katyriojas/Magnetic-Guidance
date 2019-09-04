% Determine Calibration Curve
% 9/3/19
clear all; close all; clc;

masses = [2,5,10,20,50]*10^-3; %[kg]
weights = 9.809*masses; % [N]
weights = weights*10^3; % [mN]
Fx = [18.2,45.5,91,183,458.5]; %[mN]
Fy = [18.0,46.5,91.5,183,459]; %[mN]
Fz = [18.0,45.5,92.0,184,461]; %[mN]

figure(1); grid on; hold on;
xlabel('Measured Forces (mN)');
ylabel('Calibrated Weights (mN)');
title('Calibrated Weights vs. Measured Forces');
plot(Fx,weights,'r*-','LineWidth',0.1);
plot(Fy,weights,'b*-','LineWidth',0.1);
plot(Fz,weights,'k*-','LineWidth',0.1);

Ax = Fx'\weights';
Ay = Fy'\weights';
Az = Fz'\weights';

% Mean slope
avgSlope = mean([Ax,Ay,Az]);

text(250,200,strcat('Slope Fx: ',' ',num2str(Ax)));
text(250,150,strcat('Slope Fy: ',' ',num2str(Ay)));
text(250,100,strcat('Slope Fz: ',' ',num2str(Az)));
text(250,50,strcat('Avg Slope:',' ',num2str(avgSlope)));

% Now let's plot the scaled curves
plot(Fx*avgSlope,weights,'m*--','LineWidth',0.1);
plot(Fy*avgSlope,weights,'c*--','LineWidth',0.1);
plot(Fz*avgSlope,weights,'Color',[0.7,0.75,0.71],'LineStyle','--','Marker','*','LineWidth',0.1);

legend('Fx measured','Fy measured','Fz measured', 'Fx cal', 'Fy cal','Fz cal');

