classdef SmaractData
   properties
       path
       raw
       i_start
       i_end       
   end
   properties (Dependent, SetAccess = private)
       time_unix % [ns->UNIX time] Nx1
       time_start_unix  % [ns->UNIX time] time when movement starts
       time_end_unix    % [ns->UNIX time] time when movement ends
       time_moving_unix % [ns->UNIX time] Mx1
       time_moving    % [s] Mx1
       ch0 % [mm] Nx1
       ch1 % [mm] Nx1
       ch0_moving % [mm] Mx1 ch0 positions while moving
       ch1_moving % [mm] Mx1 ch0 positions while moving
       i_range
   end
   
   methods
       function obj = SmaractData(filepath)
           obj.path = filepath;
           
           % import CSV
           obj.raw = importRosSmaractCSV(filepath);
           
           % smaract indices where movement starts/ends
           obj.i_start = find(diff(obj.ch0) > 0.001, 1);
           obj.i_end = length(obj.ch0) - find( abs(diff(flipud(obj.ch0))) > 0.001, 1);
       end
       
       function time_unix = get.time_unix(obj)
           time_unix = obj.raw.time;       
       end
       
       function time_moving_unix = get.time_moving_unix(obj)
           time_moving_unix = obj.time_unix(obj.i_range);
       end
       
       function time_start_unix = get.time_start_unix(obj)
           time_start_unix = obj.time_unix(obj.i_start);
       end

       function time_end_unix = get.time_end_unix(obj)
           time_end_unix = obj.time_unix(obj.i_end);
       end

       function time_moving = get.time_moving(obj)
           time_moving = (obj.time_unix - obj.time_start_unix)/1e9;
       end
       
       function ch0 = get.ch0(obj)
           ch0 = obj.raw.ch0;
       end
       
       function ch1 = get.ch1(obj)
           ch1 = obj.raw.ch1;
       end
       
       function ch0_moving = get.ch0_moving(obj)
          ch0_moving = obj.ch0(obj.i_range);
       end
 
       function ch1_moving = get.ch1_moving(obj)
          ch1_moving = obj.ch1(obj.i_range);
       end
       
       function i_range = get.i_range(obj)
           i_range = [obj.i_start:obj.i_end]';
       end
   end
    
end