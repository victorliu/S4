-- Duplicates Table 1, M2G2 structure, in
-- "Multilayer films composed of periodic magneto-optical and dielectric layers for use as Faraday rotators"
-- S. Sakaguchi and N. Sugimoto
-- Optics Communications 162, p. 64-70 (1999)

for n = 5,10 do
	S = S4.NewSimulation()
	S:SetLattice({1,0}, {0,1})
	S:SetNumG(1)

	-- Material definition

	S:AddMaterial("MO", {
		{4.75,0}, {0,0.00269}, {0,0},
		{0,-0.00269}, {4.75,0}, {0,0},
		{0, 0}, {0, 0}, {4.75,0}
		})

	S:AddMaterial("Dielectric", {2.5,0})
	S:AddMaterial("Vacuum", {1,0})

	S:AddLayer('AirAbove', 0, 'Vacuum')
	S:AddLayer('mlayer1', 1, 'MO')
	S:AddLayer('dlayer1', 1.3784, 'Dielectric')
	for i=2,n do
		S:AddLayerCopy('mlayer'..i, 1, 'mlayer1')
		S:AddLayerCopy('dlayer'..i, 1.3784, 'dlayer1')
	end
	S:AddLayerCopy('m2', 2, 'mlayer1')
	for i=1,2*n do
		S:AddLayerCopy('dlayer'..i+n, 1.3784, 'dlayer1')
		S:AddLayerCopy('mlayer'..i+n, 1, 'mlayer1')
	end
	S:AddLayerCopy('d2', 2.7568, 'dlayer1')
	for i=1,n do
		S:AddLayerCopy('mlayer'..i+3*n, 1, 'mlayer1')
		S:AddLayerCopy('dlayer'..i+3*n, 1.3784, 'dlayer1')
	end
	S:AddLayerCopy('AirBelow', 0, 'AirAbove') -- layer to copy

	S:SetExcitationPlanewave(
		{0,0}, -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
		{0,0}, -- s-polarization amplitude and phase (in degrees)
		{1,0}) -- p-polarization amplitude and phase

	S:SetFrequency(0.25/math.sqrt(4.75))
	t = S:GetPoyntingFlux('AirBelow', 0)
	Exr,Eyr,Ezr,Exi,Eyi,Ezi = S:GetEField({0,0,100})
	Ex = math.sqrt(Exr*Exr + Exi*Exi)
	Ey = math.sqrt(Eyr*Eyr + Eyi*Eyi)
	print(n, t, math.deg(math.atan2(Ey,Ex)))
end
