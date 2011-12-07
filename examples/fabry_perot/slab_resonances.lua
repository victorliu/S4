-- Outputs the quasinormal modes of the structure, in the form
-- of the poles of the system scattering matrix. One can work
-- out by hand where the resonances should be, and find that it
-- agrees well with the results here. We also illustrate the
-- behavior when -Im(omega) > Re(omega).

-- For this system, n = 2, and the normal incidence Fresnel
-- reflection coefficient is (n-1)/(n+1) = 1/3 = r
-- We expect resonances at
--   omega = 2pi*m/4 + 0.5*i*log(r), m is an integer
-- The frequency is scaled by 2pi internally, so in the plot we
-- should see resonances at: m/4 - i*0.0874248
-- The Q is -Re(omega_r)/(2*Im(omega_r)) where omega_r is the
-- resonance frequency.

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

for fi=0,-0.2,-0.002 do
	for fr=0.1,1.0,0.01 do
		S:SetFrequency(fr, fi)
		mantr,manti,base,expo = S:GetSMatrixDeterminant()
		-- we know for a simple structure like this we won't get overflow:
		mantr = mantr*math.pow(base,expo)
		manti = manti*math.pow(base,expo)
		manta = math.sqrt(mantr*mantr + manti*manti)
		print (fr .. '\t' .. fi .. '\t' .. mantr .. '\t' .. manti .. '\t' .. manta);
	end
	print('')
end
