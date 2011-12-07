-- In a 1D pattern, the pattern should be specified only with rectangles.
-- The half-height of the rectangle should be 0.5/L where L is the lattice
-- constant. (Assume a unit cell of unit area).

L = 1
S = S4.NewSimulation()
S:SetLattice({L,0}, {0,0}) -- 1D lattice
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
	                       {0.25, 0.5/L}) -- half-widths
S:AddLayerCopy('AirBelow', -- new layer name
               0,          -- thickness
               'AirAbove') -- layer to copy

-- E polarized along the grating periodicity direction
S:SetExcitationPlanewave(
	{0,0},  -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
	{1,0},  -- s-polarization amplitude and phase (in degrees)
	{0,0})  -- p-polarization amplitude and phase

G = S:GetGList()
io.stdout:write('#')
for i,g in ipairs(G) do
	io.stdout:write("\t(" .. g[1] .. ',' .. g[2] .. ')')
end
io.stdout:write("\n")


for freq=0.25,3.2,0.005 do
	S:SetFrequency(freq)
	
	-- backward should be zero
	power = S:GetPoyntingFluxByOrder('AirBelow', -- layer in which to get
		                         0)          -- z-offset
	io.stdout:write(freq)
	for i,v in pairs(power) do
		io.stdout:write("\t" .. v[1])
	end
	io.stdout:write("\n")
end
