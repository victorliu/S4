package.cpath = '../../../build;' .. package.cpath
S4 = require('S4v2')

-- Bottom pane of Fig. 12 in
-- Shanhui Fan and J. D. Joannopoulos,
-- "Analysis of guided resonances in photonic crystal slabs",
-- Phys. Rev. B, Vol. 65, 235112

S = S4.NewSimulation{ lattice = {{1,0}, {0,1}}, bases = 100 }
materials = {}
materials.Silicon = S:AddMaterial{ epsilon = {12,0} }
materials.Vacuum  = S:AddMaterial{ epsilon = { 1,0} }

layers = {}
layers[1] = S:AddLayer{ thickness = 0, material = materials.Vacuum }
layers[2] = S:AddLayer{ thickness = 0.5, material = materials.Silicon}
layers[2]:SetRegion{ shape = 'circle', radius = 0.2, center = {0, 0}, material = materials.Vacuum }
layers[3] = S:AddLayer{ thickness = 0, copy = layers[1] }

S:ExcitationPlanewave{ k = { 0, 0, 1 }, u = { 1, 0, 0 }, cu = { 1, 0 }, cv = { 0, 0 } }

for freq = 0.25, 0.6, 0.005 do
	S:SetFrequency(freq)
	forward, backward = layers[1]:GetPowerFlux()
	forward = layers[3]:GetPowerFlux()
	print (freq .. '\t' .. forward .. '\t' .. backward)
	io.stdout:flush()
end
