-- In a 1D pattern, the pattern should be specified only with rectangles.
-- The half-height of the rectangle should be 0.5/L where L is the lattice
-- constant. (Assume a unit cell of unit area).

pcall(loadstring(S4.arg))

if not NG then NG = 10 end

L = 1
S = S4.NewSimulation()
S:SetLattice({L,0}, {0,0}) -- 1D lattice
S:SetNumG(NG)

-- Material definition
S:AddMaterial("Metal", {-30,1}) -- real and imag parts
S:AddMaterial("Vacuum", {1,0})

S:AddLayer(
	'AirAbove', --name
	0,          --thickness
	'Vacuum')   --background material
S:AddLayer('Slab', 0.5, 'Vacuum')
S:SetLayerPatternRectangle('Slab',        -- which layer to alter
                           'Metal',     -- material in rectangle
	                       {0,0},         -- center
	                       0,             -- tilt angle (degrees)
	                       {0.25, 0.5/L}) -- half-widths
S:AddLayerCopy('AirBelow', -- new layer name
               0,          -- thickness
               'AirAbove') -- layer to copy

-- E polarized along the grating periodicity direction
if pol == 'p' then
	S:SetExcitationPlanewave(
		{0,0},  -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
		{0,0},  -- s-polarization amplitude and phase (in degrees)
		{1,0})  -- p-polarization amplitude and phase
else
	S:SetExcitationPlanewave(
		{0,0},  -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
		{1,0},  -- s-polarization amplitude and phase (in degrees)
		{0,0})  -- p-polarization amplitude and phase
end

if form == 'ref' then
	S:UseDiscretizedEpsilon(false)
elseif form == 'lan' then
	S:UseLanczosSmoothing(true)
elseif form == 'fft' then
	S:UseDiscretizedEpsilon(true)
	S:SetResolution(8)
elseif form == 'sps' then
	S:UseDiscretizedEpsilon(true)
	S:UseSubpixelSmoothing()
	S:SetResolution(8)
elseif form == 'nv' then
	S:UsePolarizationDecomposition()
	S:UseNormalVectorBasis()
	S:SetResolution(8)
elseif form == 'cpx' then
	S:UsePolarizationDecomposition()
	S:UseJonesVectorBasis()
	S:SetResolution(8)
elseif form == 'pol' then
	S:UsePolarizationDecomposition()
	S:SetResolution(8)
end

sampler = S4.NewSpectrumSampler(0.01, 2.5, {
	InitialNumPoints = 33,
	RangeThreshold = 0.001,
	MaxBend = math.cos(math.rad(10)),
	MinimumSpacing = 1e-6
	})
while not sampler:IsDone() do
	local freq = sampler:GetFrequency()
	S:SetFrequency(freq)
	-- backward should be zero
	local forward,backward = S:GetPoyntingFlux('AirBelow', -- layer in which to get
		                  0)          -- z-offset
	
	print(freq .. '\t' .. forward);

	sampler:SubmitResult(forward)
end

