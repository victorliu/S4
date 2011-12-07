-- Fig. 2a of
-- Victor Liu, Michelle Povinelli, and Shanhui Fan,
-- "Resonance-enhanced optical forces between coupled photonic crystal slabs",
-- Optics Express, Vol. 17, No. 24, 2009
-- A factor of 0.5 due to time averaging is not taken into account here


S = S4.NewSimulation()
S:SetLattice({1,0}, {0,0}) -- 1D lattice
S:SetNumG(27)

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
	                       {0.25, 0.5}) -- half-widths
S:AddLayerCopy('Spacer', 0.5, 'AirAbove')
S:AddLayer('Slab2', 0.5, 'Vacuum')
S:SetLayerPatternRectangle('Slab2',       -- which layer to alter
                           'Silicon',     -- material in rectangle
	                       {0.15,0},      -- center
	                       0,             -- tilt angle (degrees)
	                       {0.25, 0.5}) -- half-widths
S:AddLayerCopy('AirBelow', 0, 'AirAbove')

-- E polarized along the grating "rods"
S:SetExcitationPlanewave(
	{0,0},  -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
	{0,0},  -- s-polarization amplitude and phase (in degrees)
	{1,0})  -- p-polarization amplitude and phase

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
