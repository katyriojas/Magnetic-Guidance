classdef Nano17Data
   properties
       path
       raw
   end
   
   properties (Dependent, SetAccess = private)
        time  % [ns] Nx1, times corresponding to each force measurement during insertion
        force % [mN] Nx3, force measurements taken during insertion
        Fmag % [mN] Nx1, ||force_insertion||
        Fx   % [mN] Nx1, force_insertion(:,1)
        Fy   % [mN] Nx1, force_insertion(:,2)
        Fz   % [mN] Nx1, force_insertion(:,3)
        torque % [Nmm] Nx3
        Tx      % [Nmm] Nx1
        Ty
        Tz
        force_insertion_smooth
        Fmag_smooth
        Fx_smooth
        Fy_smooth
        Fz_smooth
   end
    
   properties (Access = private)
        smooth_span; % proportion of points to use for smoothing (default = 0.06)
        force_insertion_smooth_ % only recomputed if smooth_span has changed
   end
    
   methods
       function obj = Nano17Data(filepath)
           obj.path = filepath;
           
           % import CSV
           obj.raw = importRosNano17CSV(filepath);

           obj = obj.setSmoothSpan(0.06);
       end
       
       function obj = setSmoothSpan(obj, smooth_span)
            % recompute if smooth_span is changed
            if isempty(obj.smooth_span) || (smooth_span ~= obj.smooth_span)
                obj.smooth_span = smooth_span;
                
                obj.force_insertion_smooth_ = ...
                    [ smooth(obj.time, obj.Fx, smooth_span, 'loess'), ...
                      smooth(obj.time, obj.Fy, smooth_span, 'loess'), ...
                      smooth(obj.time, obj.Fz, smooth_span, 'loess')];
            end
       end
        
       function time = get.time(obj)
           time = obj.raw.time;       
       end
       
       function forces = get.force(obj)
           forces = [ obj.Fx, obj.Fy, obj.Fz];
       end
       
       function torques = get.torque(obj)
           torques = [ obj.Tx, obj.Ty, obj.Tz];
       end
       
       function Fx = get.Fx(obj)
          Fx = 1000 * obj.raw.Fx; 
       end
       
       function Fy = get.Fy(obj)
          Fy = 1000 * obj.raw.Fy; 
       end
       
       function Fz = get.Fz(obj)
          Fz = 1000 * obj.raw.Fz; 
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