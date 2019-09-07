classdef Nano17Data
   properties
       path
       raw
   end
   properties (Dependent, SetAccess = private)
       time_unix % [ns] Nx1
       time   % [s] Nx1
       time_vector % [s] Nx1 subtract off first time stamp
       force % [mN] Nx3
       Fmag   % [mN] Nx1
       Fx     % [mN] Nx1
       Fy
       Fz
       torque % [Nmm] Nx3
       Tx      % [Nmm] Nx1
       Ty
       Tz
       force_insertion_smooth
       Fx_smooth
       Fy_smooth
       Fz_smooth
       Fmag_smooth
   end
   
   properties (SetAccess = private)
       smooth_span
       cal_slope
       force_insertion_smooth_
   end
   
   methods
       function obj = Nano17Data(filepath,cal_slope)
           obj.path = filepath;
           
           % import CSV
           obj.raw = importRosNano17CSV(filepath);
           
           % set cal_slope
           obj.cal_slope = cal_slope;
           
       end
       
       function obj = setSmoothSpan(obj, smooth_span)
            % recompute if smooth_span is changed
            if isempty(obj.smooth_span) || (smooth_span ~= obj.smooth_span)
                obj.smooth_span = smooth_span;
                obj.force_insertion_smooth_ = ...
                    [smooth(obj.time, obj.Fx, obj.smooth_span, 'loess'), ...
                     smooth(obj.time, obj.Fy, obj.smooth_span, 'loess'), ...
                     smooth(obj.time, obj.Fz, obj.smooth_span, 'loess')];
            end
       end
       
       function time_unix = get.time_unix(obj)
           time_unix = obj.raw.time;
       end
       
       function time = get.time(obj)
           time = 10^-9*obj.raw.time;
       end
       
       function time_vector = get.time_vector(obj)
           time_vector = obj.time - obj.time(1); % cumulative time vec
       end
       
       function forces = get.force(obj)
           forces = obj.cal_slope * [obj.Fx, obj.Fy, obj.Fz];
       end
       
       function torques = get.torque(obj)
           torques = [obj.Tx, obj.Ty, obj.Tz];
       end
       
       function Fx = get.Fx(obj)
          Fx = obj.cal_slope * 1000 * obj.raw.Fx; 
       end
       
       function Fy = get.Fy(obj)
          Fy = obj.cal_slope * 1000 * obj.raw.Fy; 
       end
       
       function Fz = get.Fz(obj)
          Fz = obj.cal_slope * 1000 * obj.raw.Fz; 
       end
       
       function Tx = get.Tx(obj)
          Tx = obj.raw.Tx; 
       end
       
       function Ty = get.Ty(obj)
          Ty = obj.raw.Ty; 
       end
       
       function Tz = get.Tz(obj)
          Tz = obj.raw.Tz; 
       end
       
       function Fmag = get.Fmag(obj)
          Fmag = sqrt(sum(obj.force.^2, 2));
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
        
   end
    
end