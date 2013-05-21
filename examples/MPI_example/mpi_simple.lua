-- In a 1D pattern, the pattern should be specified only with rectangles.
-- The y-dimension of the rectangles is ignored.

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
	                       {0.25, 0}) -- half-widths
S:AddLayerCopy('AirBelow', -- new layer name
               0,          -- thickness
               'AirAbove') -- layer to copy

-- E polarized along the grating periodicity direction
S:SetExcitationPlanewave(
	{0,0},  -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
	{0,0},  -- s-polarization amplitude and phase (in degrees)
	{1,0})  -- p-polarization amplitude and phase

--S:UsePolarizationDecomposition()

freq_start = 0.25
freq_end = 0.7
freq_step = 0.005

my_freq_start = freq_start + S4.MPIRank/S4.MPISize * (freq_end - freq_start)
my_freq_end   = my_freq_start + (freq_end - freq_start) / S4.MPISize
outfile = io.open('out' .. S4.MPIRank .. '.txt', 'w')

for freq=my_freq_start,my_freq_end,freq_step do
	S:SetFrequency(freq)
	
	-- backward should be zero
	forward,backward = S:GetPoyntingFlux('AirBelow', -- layer in which to get
		                                 0)          -- z-offset
	
	outfile:write(freq .. '\t' .. forward .. '\n');
end
