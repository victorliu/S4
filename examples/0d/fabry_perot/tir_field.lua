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
	'Dielectric')   --background material
S:AddLayer('Below', 0, 'Vacuum')

S:SetFrequency(1);


S:SetExcitationPlanewave(
	{80,0}, -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
	{1,0},     -- s-polarization amplitude and phase (in degrees)
	{0,0})     -- p-polarization amplitude and phase


for x = -0.5,0.5, 0.01 do
	for z = -0.5,0.5, 0.01 do
		Ex,Ey,Ez = S:GetEField({x,0,z})
		print(x .. '\t' .. z.. '\t' .. Ex.. '\t' .. Ey .. '\t' .. Ez);
	end
	print('');
end
