%% Create tables of peak forces and final insertion depths

%%%%%%%%%%%
% PHANTOM %
%%%%%%%%%%%

% peak forces (from raw force measurements)
for i_trial = 1:length(data_robotic_phantom_trim)
    Fpeak.phantom.manual(i_trial) = max(data_manual_phantom(i_trial).Fmag_trimmed);
    Fpeak.phantom.nomag(i_trial)  = max(data_robotic_phantom_trim(i_trial).nomag_mea.Fmag);
    Fpeak.phantom.mag(i_trial)    = max(data_robotic_phantom_trim(i_trial).mag.Fmag);
end

% final angular insertion depths
for i_trial = 1:length(data_robotic_phantom_trim)
    AID.phantom.manual(i_trial) = max(data_manual_phantom(i_trial).interp_angdepth);
    AID.phantom.nomag(i_trial)  = max(data_robotic_phantom_trim(i_trial).nomag_mea_interp_angdepth);
    AID.phantom.mag(i_trial)    = max(data_robotic_phantom_trim(i_trial).mag_interp_angdepth);
end

method_names = {'Manual'; 'Unguided'; 'Guided'};
stat_names   = {'Fpeak'; 'Fpeak_std'; 'AID'; 'AID_std'};

phantom_table = table( [mean(Fpeak.phantom.manual); mean(Fpeak.phantom.nomag);  mean(Fpeak.phantom.mag)],...
                       [std(Fpeak.phantom.manual);  std(Fpeak.phantom.nomag);   std(Fpeak.phantom.mag)],...
                       [mean(AID.phantom.manual);   mean(AID.phantom.nomag);    mean(AID.phantom.mag)],...
                       [std(AID.phantom.manual);    std(AID.phantom.nomag);     std(AID.phantom.mag)],...
                       'RowNames', method_names, 'VariableNames', stat_names);

disp(phantom_table);

%%%%%%%%%%%
% Cadaver %
%%%%%%%%%%%

% peak forces (from raw force measurements)
for i_trial = 1:length(data_robotic_cadaver)
    Fpeak.cadaver.manual(i_trial) = max(data_manual_cadaver(i_trial).Fmag_trimmed);
    Fpeak.cadaver.nomag(i_trial)  = max(data_robotic_cadaver(i_trial).nomag_Fmag_trimmed);
    Fpeak.cadaver.mag(i_trial)    = max(data_robotic_cadaver(i_trial).mag_Fmag_trimmed);
end

% final angular insertion depths (from CIP)
AID.cadaver.manual = [498, 502, 438];
AID.cadaver.nomag  = [369, 460, 482];
AID.cadaver.mag    = [322, 462, 593];


% create table for AID/Fpeak
method_names = {'Manual'; 'Unguided'; 'Guided'};
stat_names   = {'Fpeak'; 'Fpeak_std'; 'AID'; 'AID_std'};

cadaver_table = table( [mean(Fpeak.cadaver.manual); mean(Fpeak.cadaver.nomag);  mean(Fpeak.cadaver.mag)],...
                       [std(Fpeak.cadaver.manual);  std(Fpeak.cadaver.nomag);   std(Fpeak.cadaver.mag)],...
                       [mean(AID.cadaver.manual);   mean(AID.cadaver.nomag);    mean(AID.cadaver.mag)],...
                       [std(AID.cadaver.manual);    std(AID.cadaver.nomag);     std(AID.cadaver.mag)],...
                       'RowNames', method_names, 'VariableNames', stat_names);

disp(cadaver_table);

%% Bar Graph of Final Angular Depths
if exist('hf_AID','var')
    if isvalid(hf_AID)
        close(hf_AID)
    end
end

hf_AID = figure;
hold on

fs = 10; % font size
ms = 40; % marker size

title('Final Angular Insertion Depths')

% group together by method (manual/unguided/guided)
bar_data.AID = [phantom_table.AID'; cadaver_table.AID'];

% create bar graph
h_bar_AID = bar(bar_data.AID);
h_bar_AID(1).FaceColor = [1 0 0]; % manual   -> red
h_bar_AID(2).FaceColor = [0 0 1]; % unguided -> blue
h_bar_AID(3).FaceColor = [0 1 0]; % guided   -> green


% use min/max for error bars
lows  = [min(AID.phantom.manual), min(AID.phantom.nomag), min(AID.phantom.mag);...
         min(AID.cadaver.manual), min(AID.cadaver.nomag), min(AID.cadaver.mag)];
err_neg = bar_data.AID - lows;

highs = [max(AID.phantom.manual), max(AID.phantom.nomag), max(AID.phantom.mag);...
         max(AID.cadaver.manual), max(AID.cadaver.nomag), max(AID.cadaver.mag)];
err_pos = highs - bar_data.AID;


% find center of each bar
drawnow % bars must be created before we can get the locations
for k1 = 1:size(bar_data.AID,2)
    ctr(k1,:) = bsxfun(@plus, h_bar_AID(k1).XData, h_bar_AID(k1).XOffset'); % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
end

% plot individual trials
scatter(repmat(ctr(1),[1 4]), AID.phantom.manual, ms, 'k', 'o', 'LineWidth',1);
scatter(repmat(ctr(2),[1 4]), AID.phantom.nomag,  ms, 'k', 'o', 'LineWidth',1);
scatter(repmat(ctr(3),[1 4]), AID.phantom.mag,    ms, 'k', 'o', 'LineWidth',1);
scatter(repmat(ctr(4),[1 3]), AID.cadaver.manual, ms, 'k', 'o', 'LineWidth',1);
scatter(repmat(ctr(5),[1 3]), AID.cadaver.nomag,  ms, 'k', 'o', 'LineWidth',1);
scatter(repmat(ctr(6),[1 3]), AID.cadaver.mag,    ms, 'k', 'o', 'LineWidth',1);

% plot error bars
% errorbar(ctr', bar_data.AID, err_neg, err_pos, 'k', 'linestyle', 'none');

% add labels etc.
xticks([1 2])
xticklabels({'Phantom', 'Cadaver'})
ylabel('Angular Insertion Depth (\circ)', 'FontWeight','bold', 'FontSize',fs, 'FontName','Times'); 
legend('Manual','Robotic','Robotic &\newlineMagnetic Steering', 'Location','nw', 'FontName','Times', 'FontSize',fs);
set(gca,'YGrid','on')


%% Bar Graph of Peak Forces
% if exist('hf_Fpeak','var')
%     if isvalid(hf_Fpeak)
%         close(hf_Fpeak)
%     end
% end
% 
% hf_Fpeak = figure;
% hold on
% 
% title('Peak Insertion Forces')
% 
% % group together by method (manual/unguided/guided)
% bar_data.phantom = [phantom_table.Fpeak'; cadaver_table.Fpeak'];
% 
% % create bar graph
% h_bar_Fpeak = bar(bar_data.phantom);
% 
% 
% % use min/max for error bars
% lows  = [min(Fpeak.phantom.manual), min(Fpeak.phantom.nomag), min(Fpeak.phantom.mag);...
%          min(Fpeak.cadaver.manual), min(Fpeak.cadaver.nomag), min(Fpeak.cadaver.mag)];
% err_neg = [phantom_table.Fpeak'; cadaver_table.Fpeak'] - lows;
% 
% highs = [max(Fpeak.phantom.manual), max(Fpeak.phantom.nomag), max(Fpeak.phantom.mag);...
%          max(Fpeak.cadaver.manual), max(Fpeak.cadaver.nomag), max(Fpeak.cadaver.mag)];
% err_pos = highs - [phantom_table.Fpeak'; cadaver_table.Fpeak'];
% 
% 
% % find center of each bar
% drawnow % bars must be created before we can get the locations
% for k1 = 1:size(bar_data.phantom,2)
%     ctr(k1,:) = bsxfun(@plus, h_bar_Fpeak(k1).XData, h_bar_Fpeak(k1).XOffset');       % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
% end
% 
% % plot error bars
% errorbar(ctr', bar_data.phantom, err_neg, err_pos, 'k', 'linestyle', 'none');                                 % Plot Error Bars
% 
% % add labels etc.
% xticks([1 2])
% xticklabels({'Phantom', 'Cadaver'})
% ylabel('||F|| (mN)','FontName','Times','FontWeight','bold');
% legend('Manual','Robotic','Robotic &\newlineMagnetic Steering', 'Location','best','FontName','Times');


