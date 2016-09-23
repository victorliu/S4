function run_sim(pvw)
	local n2 = 1.5
	local f = 3.3098*pvw^2 - 6.2731*pvw + 2.991
	local delta_n_eff = -1.9076*f^2 + 1.8060*f - 0.0214
	local d = 1/(pvw * 2 * delta_n_eff)
	local angle = math.deg(math.asin(0.5/pvw))

	S = S4.NewSimulation{ Lattice = 1, BasisSize = 27 }

	-- Material definition
	S:AddMaterial("FusedSilica", {n2^2,0}) -- real and imag parts
	S:AddMaterial("Vacuum", {1,0})

	S:AddLayer(
		'AirAbove', --name
		0,          --thickness
		'Vacuum')   --background material
	S:AddLayer('Grating', d, 'Vacuum')
	S:SetLayerPatternRectangle('Grating',     -- which layer to alter
							   'FusedSilica', -- material in rectangle
							   {0,0},         -- center
							   0,             -- tilt angle (degrees)
							   {0.5*f, 0}) -- half-widths
	S:AddLayer('Substrate', 0, 'FusedSilica')

	-- E polarized along the grating periodicity direction
	S:SetExcitationPlanewave(
		{angle,0}, -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
		{1,0},  -- TE-polarization amplitude and phase (in degrees)
		{0,0})  -- TM-polarization amplitude and phase

	S:SetFrequency(pvw)
	incident = S:GetPowerFlux('AirAbove')
	local P = S:GetPowerFluxByOrder('Substrate', 0)
	local de_TE = P[3][1]/incident


	-- E polarized along the grating periodicity direction
	S:SetExcitationPlanewave(
		{angle,0}, -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
		{0,0},  -- TE-polarization amplitude and phase (in degrees)
		{1,0})  -- TM-polarization amplitude and phase
	S:SetFrequency(pvw)
	incident = S:GetPowerFlux('AirAbove')
	local P = S:GetPowerFluxByOrder('Substrate', 0)
	local de_TM = P[1][1]/incident
	return de_TE, de_TM
end

print('pvw', 'TE', 'TM')
for pvw=0.501,0.9,0.001 do
	print(pvw, run_sim(pvw))
end
