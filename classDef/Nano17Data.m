classdef Nano17Data
   properties
       path
       raw
   end
   properties (Dependent, SetAccess = private)
       time   % [ns] Nx1
       force  % [mN] Nx3
       Fmag   % [mN] Nx1
       Fx     % [mN] Nx1
       Fy
       Fz
       torque % [Nmm] Nx3
       Tmag   % [Nmm] Nx1
       Tx     % [Nmm] Nx1
       Ty
       Tz
   end
   
   methods
       function obj = Nano17Data(filepath)
           obj.path = filepath;
           
           % import CSV
           obj.raw = importRosNano17CSV(filepath);
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
       
      function Tmag = get.Tmag(obj)
          Tmag = sqrt(sum(obj.torque.^2, 2));
       end
   end
    
end