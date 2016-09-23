-- Example 2 in
--  Lifeng Li,
--  "New formulation of the Fourier modal method for crossed surface-relief gratings"
--  Journal of the Optical Society of America A, Vol. 14, No. 10, p. 2758 (1997)
-- This is Fig. 7, curve OCA

S = S4.NewSimulation()
S:SetLattice({1,0}, {0.5,0.5*math.sqrt(3)})
S:SetNumG(100)
S:AddMaterial("Dielectric", {2.56,0}) -- real and imag parts
S:AddMaterial("Vacuum", {1,0})

S:AddLayer('StuffAbove', 0 , 'Vacuum')
S:AddLayer('Slab', 0.5, 'Vacuum')
S:SetLayerPatternCircle('Slab', 'Dielectric', {0,0}, 0.25)
S:AddLayer('StuffBelow', 0, 'Dielectric')

S:SetExcitationPlanewave(
	{30,30}, -- incidence angles
	{0,0}, -- s-polarization amplitude and phase (in degrees)
	{1,0}) -- p-polarization amplitude and phase

S:SetFrequency(2.0000001)

if     "original" == S4.arg then
elseif "new" == S4.arg then
		S:UsePolarizationDecomposition()
elseif "normal" == S4.arg then
		S:UsePolarizationDecomposition()
		S:UseNormalVectorBasis()
elseif "complex" == S4.arg then
		S:UsePolarizationDecomposition()
		S:UseJonesVectorBasis()
end
S:SetResolution(8)

--[[
for y = -1.5,1.5,0.02 do
	for x = -1.5,1.5,0.02 do
		print(x, y, S:GetEpsilon({x,y,0.25}))
	end
	print('')
end
]]

for ng = 41,401,40 do
	S:SetNumG(ng)
	power_inc = S:GetPoyntingFlux('StuffAbove', 0)
	G = S:GetGList()
	Gi = 1
	for i,g in ipairs(G) do
		if -1 == g[1] and -2 == g[2] then
			Gi = i
			break
		end
	end
	
	P = S:GetPoyntingFluxByOrder('StuffBelow', 0)

	print(S:GetNumG(), P[Gi][1]/power_inc, P[Gi][2]/power_inc)
	io.stdout:flush()
end

--Lk = S:GetReciprocalLattice()
--print(Lk[1][1], Lk[1][2], Lk[2][1], Lk[2][2])

