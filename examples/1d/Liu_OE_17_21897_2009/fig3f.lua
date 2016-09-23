-- Fig. 3f of
-- Victor Liu, Michelle Povinelli, and Shanhui Fan,
-- "Resonance-enhanced optical forces between coupled photonic crystal slabs",
-- Optics Express, Vol. 17, No. 24, 2009

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
S:AddLayerCopy('Spacer', 0.65, 'AirAbove')
S:AddLayerCopy('Slab2', 0.5, 'Slab')
S:AddLayerCopy('AirBelow', 0, 'AirAbove')

-- E polarized along the grating "rods"
S:SetExcitationPlanewave(
	{0,0},  -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
	{0,0},  -- s-polarization amplitude and phase (in degrees)
	{1,0})  -- p-polarization amplitude and phase

S:SetFrequency(0.6306178089044143)

for x=-0.5,3.5,0.02 do
	for z=-1,2.65,0.02 do
		Ex,Ey,Ez = S:GetEField({x,0,z})
		print(x..'\t'..z..'\t'..Ey)
	end
	print('')
end
