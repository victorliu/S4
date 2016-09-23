S = S4.NewSimulation{ Lattice = 1, BasisSize = 27 }

-- Material definition
S:AddMaterial("Silicon", {12,0}) -- real and imag parts
S:AddMaterial("Vacuum", {1,0})

S:AddLayer(
	'AirAbove', --name
	0,          --thickness
	'Vacuum')   --background material
S:AddLayer('Slab', 0.5, 'Vacuum')
S:SetLayerPatternRectangle('Slab',        -- which layer to alter
                           'Silicon',     -- material in rectangle
	                       {0,0},         -- center
	                       0,             -- tilt angle (degrees)
	                       {0.25, 0}) -- half-widths
S:AddLayerCopy('AirBelow', -- new layer name
               0,          -- thickness
               'AirAbove') -- layer to copy

-- E polarized along the grating periodicity direction
S:SetExcitationPlanewave(
	{0,0},  -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
	{1,0},  -- s-polarization amplitude and phase (in degrees)
	{0,0})  -- p-polarization amplitude and phase

for freq=0.25,0.7,0.005 do
	S:SetFrequency(freq)

	-- backward should be zero
	forward,backward = S:GetPowerFlux('AirBelow', -- layer in which to get
		                                 0)          -- z-offset

	print(freq .. '\t' .. forward);
end
