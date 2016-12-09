classdef S4Simulation < handle
	properties (SetAccess = private, Hidden = true)
		S;
	end
	methods
		function obj = S4Simulation(Lattice, Bases)
			if nargin > 0
				obj.S = S4mex('new', Lattice, Bases);
				S4Material();
				S4Layer();
			end
		end
		function delete(obj)
			S4mex('destroy', obj.S);
		end
		function M = setMaterial(obj, epsilon)
			m = S4mex('set_material', obj.S, epsilon);
			M = S4Material(obj.S, m);
		end
		function L = addLayer(obj, thickness, material)
			l = S4mex('add_layer', obj.S, thickness, material.M);
			L = S4Layer(obj.S, l);
		end
		function setPlanewave(obj, k, u, cu, cv)
			S4mex('set_planewave', obj.S, k, u, cu, cv);
		end
		function setFrequency(obj, freq)
			S4mex('set_frequency', obj.S, freq);
		end
	end
end
