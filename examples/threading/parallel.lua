Sa = S4.NewSimulation()
Sa:SetLattice({1,0}, {0,1})
Sa:SetNumG(41)
Sa:AddMaterial("Silicon", {12,0}) -- real and imag parts
Sa:AddMaterial("Vacuum", {1,0})

Sa:AddLayer('AirAbove', 0 , 'Vacuum')
Sa:AddLayer('Slab', 0.5, 'Silicon')
Sa:SetLayerPatternCircle('Slab', 'Vacuum', {0,0}, 0.2)
Sa:AddLayerCopy('AirBelow', 0, 'AirAbove')

Sa:SetExcitationPlanewave(
	{0,0}, -- incidence angles
	{1,0}, -- s-polarization amplitude and phase (in degrees)
	{0,0}) -- p-polarization amplitude and phase

Sa:UsePolarizationDecomposition()
Sa:SetResolution(8)

Sb = Sa:Clone()
Sc = Sa:Clone()
Sd = Sa:Clone()

delta = 0.01
for freq=0.25,0.6,delta do
	fa = freq
	fb = fa + 0.25*delta
	fc = fb + 0.25*delta
	fd = fc + 0.25*delta
	Sa:SetFrequency(fa)
	Sb:SetFrequency(fb)
	Sc:SetFrequency(fc)
	Sd:SetFrequency(fd)
	
	S4.SolveInParallel('AirBelow', Sa, Sb, Sc, Sd);

	print(fa, Sa:GetPoyntingFlux('AirBelow', 0))
	print(fb, Sb:GetPoyntingFlux('AirBelow', 0))
	print(fc, Sc:GetPoyntingFlux('AirBelow', 0))
	print(fd, Sd:GetPoyntingFlux('AirBelow', 0))
end
