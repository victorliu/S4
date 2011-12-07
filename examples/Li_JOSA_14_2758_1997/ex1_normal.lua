-- Extension of Example 1A in
--  Lifeng Li,
--  "New formulation of the Fourier modal method for crossed surface-relief gratings"
--  Journal of the Optical Society of America A, Vol. 14, No. 10, p. 2758 (1997)
-- This would be in Fig. 6

S = S4.NewSimulation()
S:SetLattice({2.5,0}, {0,2.5})
S:SetNumG(200)
S:AddMaterial("Dielectric", {2.25,0}) -- real and imag parts
S:AddMaterial("Vacuum", {1,0})

S:AddLayer('StuffAbove', 0 , 'Dielectric')
S:AddLayer('Slab', 1.0, 'Vacuum')
S:SetLayerPatternRectangle('Slab', 'Dielectric', {-1.25/2,-1.25/2}, 0, {1.25/2, 1.25/2})
S:SetLayerPatternRectangle('Slab', 'Dielectric', { 1.25/2, 1.25/2}, 0, {1.25/2, 1.25/2})
S:AddLayer('AirBelow', 0, 'Vacuum')

S:SetExcitationPlanewave(
	{0,0}, -- incidence angles
	{1,0}, -- s-polarization amplitude and phase (in degrees)
	{0,0}) -- p-polarization amplitude and phase

S:UsePolarizationDecomposition()
S:UseNormalVectorBasis()
S:SetResolution(8)

--S:OutputLayerPatternRealization('Slab', 128, 128)

S:SetFrequency(1)

for ng = 41,401,40 do
	S:SetNumG(ng)
	power_inc = S:GetPoyntingFlux('StuffAbove', 0)
	G = S:GetGList()
	P = S:GetPoyntingFluxByOrder('AirBelow', 0)
	actualg = S:GetNumG()
	print(actualg, P[6][1]/power_inc) -- 1+4+1
end
