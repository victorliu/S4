
S = S4.NewSimulation()
S:SetLattice({1,0}, {0,1})
S:SetNumG(1)

-- Material definition
S:AddMaterial("Dielectric", {4,0}) -- real and imag parts
S:AddMaterial("Vacuum", {1,0})

S:AddLayer(
	'AirAbove', --name
	0,          --thickness
	'Vacuum')   --background material

S:AddLayer('Slab', 1, 'Dielectric')
S:AddLayerCopy('AirBelow', -- new layer name
               0,          -- thickness
               'AirAbove') -- layer to copy

S:SetExcitationPlanewave(
	{0,0}, -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
	{0,0}, -- s-polarization amplitude and phase (in degrees)
	{1,0}) -- p-polarization amplitude and phase

for freq=0.001,1.0,0.005 do
	S:SetFrequency(freq)
	
	-- backward should be zero
	forward,backward = S:GetPoyntingFlux('AirBelow', -- layer in which to get
		                                 0)          -- z-offset
	slab_forward,slab_backward = S:GetPoyntingFlux('Slab', 0)
	E2 = S:GetLayerElectricEnergyDensityIntegral('Slab');
	
	print(freq .. '\t' .. forward .. '\t' .. slab_forward .. '\t' .. slab_backward .. '\t' .. E2);
end

-- The energy density ought to match the analytic expression:
-- Plot[
--  n^2 (2/(n + 1))^2/(
--    1 + ((n - 1)/(n + 1))^4 - 
--     2 ((n - 1)/(n + 1))^2 Cos[
--       2 n (2 \[Pi] k)]) ((1 + ((n - 1)/(n + 1))^2) + 
--      2 Sinc[n (2 \[Pi] k)] (n - 1)/(n + 1)
--        Cos[n (2 \[Pi] k)]) /. {n -> 2}
--  , {k, 0, 1}]
