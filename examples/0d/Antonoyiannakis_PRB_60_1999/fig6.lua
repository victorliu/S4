-- Reproduces one of the curves in Fig. 6 in
-- M. I. Antonoyiannakis and J. B. Pendry,
-- "Electromagnetic forces in photonic crystals",
-- Phys. Rev. B, Vol. 60, No. 4, 1999.

S = S4.NewSimulation()
S:SetLattice({1,0}, {0,1})
S:SetNumG(1)
S:AddMaterial('Eps1', {10,0})
S:AddMaterial('Eps2', {1,0})
S:AddMaterial('Eps3', {9,0})

S:AddLayer('Above', 0 , 'Eps1')
S:AddLayer('Slab', 1, 'Eps2')
S:AddLayer('Below', 0, 'Eps3')

S:SetFrequency(0.01) -- lambda/d = 100

for phi=0,89,0.5 do
	S:SetExcitationPlanewave(
		{phi,0}, -- incidence angles
		{1,0},   -- s-polarization amplitude and phase (in degrees)
		{0,0})   -- p-polarization amplitude and phase
	
	T1x,T1y,T1z = S:GetStressTensorIntegral('Slab', 0.9)
	
	-- We use -T1z since the normal vector to the enclosing surface is -z
	print (math.cos(math.rad(phi)) .. '\t' .. -T1z)
end

