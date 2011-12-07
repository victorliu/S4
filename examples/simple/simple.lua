-- Almost the simplest example: a simple planewave travelling through vacuum

S = S4.NewSimulation()
S:SetLattice({1,0}, {0,1}) -- 1D lattice
S:SetNumG(1)

S:AddMaterial("Vacuum", {1,0})

S:AddLayer(
	'Front', --name
	0,          --thickness
	'Vacuum')   --background material
S:AddLayerCopy('Back',  -- new layer name
               0,       -- thickness
               'Front') -- layer to copy

S:SetExcitationPlanewave(
	{0,0},  -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
	{1,0},  -- s-polarization amplitude and phase (in degrees)
	{0,0})  -- p-polarization amplitude and phase

if S4.arg == nil then S:SetFrequency(0.5)
else S:SetFrequency(S4.arg)
end

for z=-2,2,0.1 do
	print(S:GetEField({0,0,z}))
	print(S:GetPoyntingFlux('Front',z))
	print(S:GetEField({0,0,z}))
	print(S:GetPoyntingFlux('Back',z))
end
