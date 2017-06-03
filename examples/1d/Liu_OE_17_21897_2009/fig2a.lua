-- Fig. 2a of
-- Victor Liu, Michelle Povinelli, and Shanhui Fan,
-- "Resonance-enhanced optical forces between coupled photonic crystal slabs",
-- Optics Express, Vol. 17, No. 24, 2009
-- A factor of 0.5 due to time averaging is not taken into account here

S4 = require('S4v2')

S = S4.NewSimulation{ lattice = 1, bases = 27 }

-- Material definition
materials = {}
materials.Vacuum  = S:AddMaterial{ epsilon = 1  }
materials.Silicon = S:AddMaterial{ epsilon = 12 }

layers = {}
layers[1] = S:AddLayer{ thickness = 0, material = materials.Vacuum }
layers[2] = S:AddLayer{ thickness = 0.5, material = materials.Vacuum }
layers[2]:SetRegion{ material = materials.Silicon, center = 0, halfwidths = 0.25 }
layers[3] = S:AddLayer{ thickness = 0.5, material = materials.Vacuum }
layers[4] = S:AddLayer{ thickness = 0.5, material = materials.Vacuum }
layers[4]:SetRegion{ material = materials.Silicon, center = 0.15, halfwidths = 0.25 }
layers[4] = S:AddLayer{ thickness = 0, material = materials.Vacuum }

S:ExcitationPlanewave{ k = {0,0,1}, u = {1,0,0}, cu = {0,0}, cv = {1,0} }

S:SetFrequency(0.57)

for dx=-0.5,0.5,0.01 do
	S:SetLayer('Slab2', 0.5)
	S:SetLayerPatternRectangle('Slab2',
                           'Silicon',     -- material in rectangle
	                       {dx,0},      -- center
	                       0,             -- tilt angle (degrees)
	                       {0.25, 0.5}) -- half-widths
	T1x,T1y,T1z = S:GetStressTensorIntegral('Spacer', 0.5)
	T2x,T2y,T2z = S:GetStressTensorIntegral('AirBelow', 0)

	print(dx..'\t'..T2x-T1x..'\t'..T2z-T1z)
end
