S4 = require('S4v2')

S = S4.NewSimulation{ lattice = 1, bases = 27 }

materials = {}
materials.Silicon = S:AddMaterial{ epsilon = 12 }
materials.Vacuum  = S:AddMaterial{ epsilon = 1  }

layers = {}
layers[1] = S:AddLayer{ thickness = 0, material = materials.Vacuum }
layers[2] = S:AddLayer{ thickness = 0.5, material = materials.Vacuum }
layers[2]:SetRegion{ material = materials.Silicon, center = 0, halfwidths = {0.25,0} }
layers[3] = S:AddLayer{ thickness = 0, material = materials.Vacuum }

S:ExcitationPlanewave{ k = {0,0,1}, u = {1,0,0}, cu = {0,0}, cv = {1,0} }

for freq=0.25,0.7,0.005 do
	S:SetFrequency(freq)

	-- backward should be zero
	inc = layers[1]:GetPowerFlux()
	forward,backward = layers[3]:GetPowerFlux()
	print(freq .. '\t' .. forward/inc);
end
