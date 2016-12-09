classdef S4Material < handle
	properties (GetAccess = {?S4Simulation, ?S4Layer})
		S;
		M;
	end
	methods
		function obj = S4Material(s, m)
			if nargin > 0
				obj.S = s;
				obj.M = m;
			end
		end
	end
end
