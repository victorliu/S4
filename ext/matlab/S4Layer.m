classdef S4Layer < handle
	properties (GetAccess = ?S4Simulation)
		S;
		L;
	end
	methods
		function obj = S4Layer(s, l)
			if nargin > 0
				obj.S = s;
				obj.L = l;
			end
		end
		function setRegion(obj, material, varargin)
			S4mex('set_layer_region', obj.S, obj.L, material.M, varargin{:});
		end
		function [f,b] = getPowerFlux(obj)
			[f,b] = S4mex('get_power_flux', obj.S, obj.L);
		end
	end
end
