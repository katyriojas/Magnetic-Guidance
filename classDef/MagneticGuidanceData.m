%% MagneticGuidanceData Class
%
%   Imports, stores, and processes data from magnetic guidance experiments
%
%
%   How to Use:
%   -> create CSV files for smaract and force topics from ROS bag file
%     -> use bash script 'bag2csv.sh'
%       caos@control:~/MagSteering$ ./bag2csv.sh
%
%   -> create a struct with filepaths to the CSV files
%     -> MUST use these field names: smaract, forces, (optional) dac
%       >> filepaths.smaract = 'C:\path\to\data\smaract_data_from_bag.csv';
%       >> filepaths.force  = 'C:\path\to\data\force_data_from_bag.csv';
% 
%   -> construct instance of class
%       >> trial2_nomag = MagneticGuidanceData(filepaths);
% 
%   -> set smoothing span to use for force data
%       >> trial2_nomag.setSmoothSpan(0.05);
% 
%   -> plot force vs insertion depth
%       >> plot(trial2_nomag.depth_insertion, trial2_nomag.Fmag, '--b')
%       >> plot(trial2_nomag.depth_insertion, trial2_nomag.Fmag_smooth, 'b')
%
%   Trevor Bruns
%   July 2019

%%
classdef MagneticGuidanceData
    properties
        smooth_span; % # of points to use for smoothing (default = 40)
    end
    properties (SetAccess = private)
        depth_insertion % [mm] Nx1, interpolated depth at each force measurement during insertion        
        smaract % SmaractData object
        nano    % Nano17Data object
%         dac
        force_i_start
        force_i_end
    end
    properties (Dependent, SetAccess = private)
        time_insertion_unix  % [ns->UNIX time] Nx1, times corresponding to each force measurement during insertion
        time_insertion  % [s] Nx1, times corresponding to each force measurement during insertion
        force_insertion % [mN] Nx3, force measurements taken during insertion
        Fmag % [mN] Nx1, ||force_insertion||
        Fx   % [mN] Nx1, force_insertion(:,1)
        Fy   % [mN] Nx1, force_insertion(:,2)
        Fz   % [mN] Nx1, force_insertion(:,3)
        Fmag_smooth
        Fx_smooth
        Fy_smooth
        Fz_smooth
        force_insertion_smooth % only recomputed if smooth_span has changed
        force_i_insertion % indices of force measurements taken during insertion
        torque_insertion % [mN] Nx3, torque measurements taken during insertion
        Tmag % [mN] Nx1, ||torque_insertion||
        Tx   % [mN] Nx1, torque_insertion(:,1)
        Ty   % [mN] Nx1, torque_insertion(:,2)
        Tz   % [mN] Nx1, torque_insertion(:,3)
        torque_insertion_smooth
        Tmag_smooth
        Tx_smooth
        Ty_smooth
        Tz_smooth
    end
    
    properties (Access = private)
        force_insertion_smooth_  % only recomputed if smooth_span has changed        
        torque_insertion_smooth_ % only recomputed if smooth_span has changed
    end

    methods
        function obj = MagneticGuidanceData(filepaths, cal_slopes)
            % import forces CSV
            if isfield(filepaths, 'force')
               obj.nano = Nano17Data(filepaths.force, cal_slopes);
            else
                error('''forces'' struct field not found')
            end
            
            % import smaract CSV
            if isfield(filepaths, 'smaract')
               obj.smaract = SmaractData(filepaths.smaract);
            else
                error('''smaract'' struct field not found')
            end
            
            % import dac voltages CSV
%             if isfield(filepaths, 'smaract')
%                obj.dac = importRosDacCSV(filepaths.dac);
%             else
%                 disp('''dac'' struct field not found')
%             end    

            % force indices corresponding to smaract start/stop (i.e. during insertion)
            obj.force_i_start = find(obj.nano.time_unix >= obj.smaract.time_start_unix, 1);
            obj.force_i_end   = find(obj.nano.time_unix >= obj.smaract.time_end_unix,   1);
            
            % interpolate to find ch0 position at each force measurement during insertion
            obj.depth_insertion = interp1(obj.smaract.time_unix, obj.smaract.ch0, obj.time_insertion_unix);
            
            % initialize default smoothing
            obj.smooth_span = 40; % [# samples]
        end
        
        function obj = set.smooth_span(obj, new_smooth_span)
            % recompute if smooth_span is changed
            if isempty(obj.smooth_span) || (new_smooth_span ~= obj.smooth_span)
                obj.smooth_span = new_smooth_span;
                obj.force_insertion_smooth_ = ...
                    [ smooth(obj.depth_insertion, obj.Fx, obj.smooth_span, 'loess'), ...
                      smooth(obj.depth_insertion, obj.Fy, obj.smooth_span, 'loess'), ...
                      smooth(obj.depth_insertion, obj.Fz, obj.smooth_span, 'loess')];
                  
                obj.torque_insertion_smooth_ = ...
                    [ smooth(obj.depth_insertion, obj.Tx, obj.smooth_span, 'loess'), ...
                      smooth(obj.depth_insertion, obj.Ty, obj.smooth_span, 'loess'), ...
                      smooth(obj.depth_insertion, obj.Tz, obj.smooth_span, 'loess')];
            end
        end
               
        function time_insertion_unix = get.time_insertion_unix(obj)
            time_insertion_unix = obj.nano.time_unix(obj.force_i_insertion);
        end
        
        function time_insertion = get.time_insertion(obj)
            time_insertion = (obj.time_insertion_unix - obj.time_insertion_unix(1))/1e9;
        end    
        
        function force_insertion = get.force_insertion(obj)
            force_insertion = obj.nano.force(obj.force_i_insertion);
        end
        
        function Fx = get.Fx(obj)
            Fx = obj.nano.Fx(obj.force_i_insertion);
        end
       
        function Fy = get.Fy(obj)
            Fy = obj.nano.Fy(obj.force_i_insertion); 
        end
       
        function Fz = get.Fz(obj)
            Fz = obj.nano.Fz(obj.force_i_insertion); 
        end
       
        function Fmag = get.Fmag(obj)
            Fmag = obj.nano.Fmag(obj.force_i_insertion);
        end
        
        function force_insertion_smooth = get.force_insertion_smooth(obj)
            force_insertion_smooth = obj.force_insertion_smooth_;
        end
        
        function Fx_smooth = get.Fx_smooth(obj)
            Fx_smooth = obj.force_insertion_smooth(:,1);
        end
       
        function Fy_smooth = get.Fy_smooth(obj)
            Fy_smooth = obj.force_insertion_smooth(:,2); 
        end
       
        function Fz_smooth = get.Fz_smooth(obj)
            Fz_smooth = obj.force_insertion_smooth(:,3); 
        end
       
        function Fmag_smooth = get.Fmag_smooth(obj)
            Fmag_smooth = sqrt(sum(obj.force_insertion_smooth.^2, 2));
        end
        
        function force_i_insertion = get.force_i_insertion(obj)
            force_i_insertion = [obj.force_i_start:obj.force_i_end]';
        end
        
        function torque_insertion = get.torque_insertion(obj)
            torque_insertion = obj.nano.torque(obj.force_i_insertion);
        end
        
        function Tx = get.Tx(obj)
            Tx = obj.nano.Tx(obj.force_i_insertion);
        end
       
        function Ty = get.Ty(obj)
            Ty = obj.nano.Ty(obj.force_i_insertion); 
        end
       
        function Tz = get.Tz(obj)
            Tz = obj.nano.Tz(obj.force_i_insertion); 
        end
       
        function Tmag = get.Tmag(obj)
            Tmag = obj.nano.Tmag(obj.force_i_insertion);
        end
        
        function torque_insertion_smooth = get.torque_insertion_smooth(obj)
            torque_insertion_smooth = obj.torque_insertion_smooth_;
        end
        
        function Tx_smooth = get.Tx_smooth(obj)
            Tx_smooth = obj.torque_insertion_smooth(:,1);
        end
       
        function Ty_smooth = get.Ty_smooth(obj)
            Ty_smooth = obj.torque_insertion_smooth(:,2); 
        end
       
        function Tz_smooth = get.Tz_smooth(obj)
            Tz_smooth = obj.torque_insertion_smooth(:,3); 
        end
       
        function Tmag_smooth = get.Tmag_smooth(obj)
            Tmag_smooth = sqrt(sum(obj.torque_insertion_smooth.^2, 2));
        end
        
    end
end