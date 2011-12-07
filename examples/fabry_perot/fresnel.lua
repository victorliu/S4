-- Outputs the Fresnel coefficients for different angles and polarizations
-- Note that these are the power coefficients. The plots should match
-- those on Wikipedia.

S = S4.NewSimulation()
S:SetLattice({1,0}, {0,1})
S:SetNumG(1)

-- Material definition
S:AddMaterial("Dielectric", {4,0}) -- real and imag parts
S:AddMaterial("Vacuum", {1,0})

S:AddLayer(
	'Above', --name
	0,          --thickness
	'Vacuum')   --background material
S:AddLayer('Below', 0, 'Dielectric')

S:SetFrequency(1);

for angle=0,90,1 do
	S:SetExcitationPlanewave(
		{angle,0}, -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
		{0,0},     -- s-polarization amplitude and phase (in degrees)
		{1,0})     -- p-polarization amplitude and phase
	inc1,backward1 = S:GetPoyntingFlux('Above', 0)
	forward1       = S:GetPoyntingFlux('Below', 0)
	backward1 = -backward1/inc1
	forward1  =  forward1/inc1
	
	S:SetExcitationPlanewave(
		{angle,0}, -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
		{1,0},     -- s-polarization amplitude and phase (in degrees)
		{0,0})     -- p-polarization amplitude and phase
	inc2,backward2 = S:GetPoyntingFlux('Above', 0)
	forward2       = S:GetPoyntingFlux('Below', 0)
	backward2 = -backward2/inc2
	forward2  =  forward2/inc2
	
	print(angle .. '\t' .. forward1 .. '\t' .. backward1 .. '\t' .. forward2 .. '\t' .. backward2);
end
