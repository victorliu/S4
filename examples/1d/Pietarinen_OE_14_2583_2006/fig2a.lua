function run_sim(d, x1, x2, x3, z1, z2, z3, angle)
	S = S4.NewSimulation{ Lattice = d, BasisSize = 21 }

	-- Material definition
	S:AddMaterial("SiO2", {1.5^2,0}) -- real and imag parts
	S:AddMaterial("Vacuum", {1,0})

	S:AddLayer(
		'Substrate', --name
		0,          --thickness
		'SiO2')   --background material
	S:AddLayer('Level1', z1, 'Vacuum')
	S:SetLayerPatternRectangle('Level1', 'SiO2', {0.5*(1+x1)*d,0}, 0, {0.5*(1-x1)*d, 0})
	S:AddLayer('Level2', z2-z1, 'Vacuum')
	S:SetLayerPatternRectangle('Level2', 'SiO2', {0.5*(1+x2)*d,0}, 0, {0.5*(1-x2)*d, 0})
	S:AddLayer('Level3', z3-z2, 'Vacuum')
	S:SetLayerPatternRectangle('Level3', 'SiO2', {0.5*(1+x3)*d,0}, 0, {0.5*(1-x3)*d, 0})
	S:AddLayer('Air', 0, 'Vacuum')

	-- E polarized along the grating periodicity direction
	S:SetExcitationPlanewave(
		{angle,0},  -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
		{1,0},  -- s-polarization amplitude and phase (in degrees)
		{0,0})  -- p-polarization amplitude and phase

	for wl=1.001,2,0.05 do
		S:SetFrequency(1/wl)

		local incident = S:GetPowerFlux('Substrate')

		local P = S:GetPowerFluxByOrder('Air', 0)
		local de = P[2][1]/incident
		print(wl, de)
	end
end

--- Table 1
run_sim(11, 0.250, 0.500, 0.750, 0.722, 1.444, 2.166, 0) -- regular profile
run_sim(11, 0.317, 0.601, 0.914, 0.784, 1.637, 2.563, 0) -- maximum efficiency
run_sim(11, 0.356, 0.610, 0.909, 0.714, 1.507, 3.815, -9.9) -- flattened efficiency
run_sim( 5, 0.250, 0.500, 0.750, 0.722, 1.444, 2.166, 0) -- regular profile
run_sim( 5, 0.337, 0.558, 0.897, 0.796, 1.200, 3.003, 1.3) -- maximum efficiency
run_sim( 5, 0.031, 0.445, 0.940, 0.203, 1.055, 6.932, 0) -- flattened efficiency
