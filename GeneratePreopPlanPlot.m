% Load Workspace data for cadaver trials preop plan
% Plot the magnet field and currents for just the magnetic steering region
%Katy Riojas
% 9/10/19
% 

clear all; close all;
load('preoperative plans\PreopPlan_2019-08-23_15-12\workspace.mat');
%%
figure(1); clf(1);
subplot(2,1,1);
grid on; hold on;
xvec_zeros = linspace(0,insertion_depth(1),101);
I_zeros = zeros(1,100);
Ix_zeros = [I_zeros,currents_scaled(1,planned_pnts(1))];
Iy_zeros = [I_zeros,currents_scaled(2,planned_pnts(1))];
Iz_zeros = [I_zeros,currents_scaled(3,planned_pnts(1))];

h1 = plot(insertion_depth,currents_scaled(1,planned_pnts),'LineWidth',1,'Color','b');
plot(xvec_zeros,Ix_zeros,'LineWidth',1,'Color','b');
h2 = plot(insertion_depth,currents_scaled(2,planned_pnts),'LineWidth',1,'Color','y');
plot(xvec_zeros,Iy_zeros,'LineWidth',1,'Color','y');
h3 = plot(insertion_depth,currents_scaled(3,planned_pnts),'LineWidth',1,'Color','r');
plot(xvec_zeros,Iz_zeros,'LineWidth',1,'Color','r');

xlabel('Insertion Depth (mm)'); ylabel('Planned Currents (A)');
legend([h1,h2,h3],'Ix','Iy','Iz');
xlim([0,27]);
xticks([0,3,6,9,12,15,18,21,24,27]);
yticks([-50,-25,0,25,50]);

% Now plot Bmag field generated
subplot(2,1,2);
grid on; hold on;
ylim([0,100]); xlim([0,27]);
yticks([0,25,50,75,100]);
xticks([0,3,6,9,12,15,18,21,24,27]);
xlabel('Insertion Depth (mm)');
ylabel('Magnetic Field Magnitude (mT)');
Bmag_zeros = [zeros(1,100),Bmag_ramp(planned_pnts(1))*1000];

plot(insertion_depth,Bmag_ramp(planned_pnts)*1000,'k','LineWidth',1);
plot(xvec_zeros,Bmag_zeros,'k','LineWidth',1);
