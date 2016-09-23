function run_angle(angle)
	S = S4.NewSimulation{ Lattice = 1.866, BasisSize = 27 }

	-- Material definition
	S:AddMaterial("FusedSilica", {1.46847^2,0}) -- real and imag parts
	S:AddMaterial("Vacuum", {1,0})

	S:AddLayer(
		'AirAbove', --name
		0,          --thickness
		'Vacuum')   --background material
	S:AddLayer('Level1', 0.333, 'Vacuum')
	S:SetLayerPatternRectangle('Level1', 'FusedSilica', {0.5*0.519,0}, 0, {0.5*0.519, 0})
	S:AddLayer('Level2', 0.281, 'Vacuum')
	S:SetLayerPatternRectangle('Level2', 'FusedSilica', {0.5*1.202,0}, 0, {0.5*1.202, 0})
	S:AddLayer('Substrate', 0, 'FusedSilica')

	-- E polarized along the grating periodicity direction
	S:SetExcitationPlanewave(
		{-angle,0},  -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
		{1,0},  -- s-polarization amplitude and phase (in degrees)
		{0,0})  -- p-polarization amplitude and phase

	S:SetFrequency(1/0.416)

	local incident = S:GetPowerFlux('AirAbove', 0)
	local P = S:GetPowerFluxByOrder('Substrate', 0)
	local T_TE = P[3][1]/incident

	-- E polarized along the grating periodicity direction
	S:SetExcitationPlanewave(
		{-angle,0},  -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
		{0,0},  -- s-polarization amplitude and phase (in degrees)
		{1,0})  -- p-polarization amplitude and phase

	S:SetFrequency(1/0.416)

	incident = S:GetPowerFlux('AirAbove', 0)
	P = S:GetPowerFluxByOrder('Substrate', 0)
	T_TM = P[3][1]/incident

	print(angle, T_TE, T_TM);
end

for angle = -60, 60, 1 do
	run_angle(angle)
end
