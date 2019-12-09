% Katy Riojas
% Process Coil Temperatures
% 12/9/19
clear all; clc;
%% Phantom Coil Temps for Guided Trials
% Folder spec
folder_path_phantom = 'data\phantom\Coil Temps\';

for ii = 1:4
    phantom_coil_temps(ii) = importCoilTemps(strcat(folder_path_phantom,'phantom_g_mea1_trial',num2str(ii),'_1.25_coiltemps.csv'));
end


%% Cadaver Coil Temps for Guided Trials
folder_path_cadaver = 'data\cadaver\Coil Temps\';
for ii = 1:3
    cadaver_coil_temps(ii) = importCoilTemps(strcat(folder_path_cadaver,'cadaver_g_mea2_trial',num2str(ii),'_coiltemps.csv'));
end

%% Plotting and Maximums
% Make two matrices where each row is a trial and each column is the max Tx,
% Ty, and Tz
max_coil_temps_phantom = zeros(4,3);
max_coil_temps_cadaver = zeros(3,3);

figure(1); clf(1); grid on; hold on;
title('Coil Temps for Phantom Trials');
xlabel('Time (s)');
ylabel('\DeltaT_{coil} (^\circC)');

for ii = 1:4
    [max_coil_temps_phantom(ii,1),max_Tx_idx] = max(phantom_coil_temps(ii).Tx);
    [max_coil_temps_phantom(ii,2),max_Ty_idx] = max(phantom_coil_temps(ii).Ty);
    [max_coil_temps_phantom(ii,3),max_Tz_idx] = max(phantom_coil_temps(ii).Tz);
   
    scatter(phantom_coil_temps(ii).time(max_Tx_idx),max_coil_temps_phantom(ii,1),'filled','r');
    scatter(phantom_coil_temps(ii).time(max_Ty_idx),max_coil_temps_phantom(ii,2),'filled','g');
    scatter(phantom_coil_temps(ii).time(max_Tz_idx),max_coil_temps_phantom(ii,3),'filled','b');
    
    plot(phantom_coil_temps(ii).time,phantom_coil_temps(ii).Tx,'r','Marker','o');
    plot(phantom_coil_temps(ii).time,phantom_coil_temps(ii).Ty,'g','Marker','o');
    plot(phantom_coil_temps(ii).time,phantom_coil_temps(ii).Tz,'b','Marker','o');
end
legend('X Coil','Y Coil','Z Coil');

figure(2); clf(2); grid on; hold on;
title('Coil Temps for Cadaver Trials');
xlabel('Time (s)');
ylabel('\DeltaT_{coil} (^\circC)');
for ii = 1:3
    [max_coil_temps_cadaver(ii,1),max_Tx_idx] = max(cadaver_coil_temps(ii).Tx);
    [max_coil_temps_cadaver(ii,2),max_Ty_idx] = max(cadaver_coil_temps(ii).Ty);
    [max_coil_temps_cadaver(ii,3),max_Tz_idx] = max(cadaver_coil_temps(ii).Tz);
    
    scatter(cadaver_coil_temps(ii).time(max_Tx_idx),max_coil_temps_cadaver(ii,1),'filled','r');
    scatter(cadaver_coil_temps(ii).time(max_Ty_idx),max_coil_temps_cadaver(ii,2),'filled','g');
    scatter(cadaver_coil_temps(ii).time(max_Tz_idx),max_coil_temps_cadaver(ii,3),'filled','b');
    
    plot(cadaver_coil_temps(ii).time,cadaver_coil_temps(ii).Tx,'r','Marker','o');
    plot(cadaver_coil_temps(ii).time,cadaver_coil_temps(ii).Ty,'g','Marker','o');
    plot(cadaver_coil_temps(ii).time,cadaver_coil_temps(ii).Tz,'b','Marker','o');
end
legend('X Coil','Y Coil','Z Coil');


%% Calculate and print maximums from each matrix
max_phantom_xyz_temps = max(max_coil_temps_phantom)
max_cadaver_xyz_temps = max(max_coil_temps_cadaver)
