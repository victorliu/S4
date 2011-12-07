pcall(loadstring(S4.arg))

if not NG then NG = 100 end

local S = S4.NewSimulation()
S:SetLattice({2,0}, {1,math.sqrt(3)})
S:SetNumG(NG)
S:AddMaterial("Dielectric", {2.56,0}) -- real and imag parts
S:AddMaterial("Vacuum", {1,0})

S:AddLayer('Before', 0 , 'Vacuum')
S:AddLayer('Slab', 1, 'Vacuum')
S:SetLayerPatternCircle('Slab', 'Dielectric', {0,0}, 0.5)
S:AddLayer('After', 0 , 'Dielectric')
S:SetFrequency(1)

S:SetExcitationPlanewave(
	{30,30}, -- incidence angles
	{0,0}, -- s-polarization amplitude and phase (in degrees)
	{1,0}) -- p-polarization amplitude and phase

if form == 'ref' then
	S:UseDiscretizedEpsilon(false)
if form == 'fft' then
	S:UseDiscretizedEpsilon(true)
	S:SetResolution(8)
elseif form == 'sps' then
	S:UseDiscretizedEpsilon(true)
	S:UseSubpixelSmoothing()
	S:SetResolution(8)
elseif form == 'nv' then
	S:UsePolarizationDecomposition()
	S:UseNormalVectorBasis()
	S:SetResolution(8)
elseif form == 'cpx' then
	S:UsePolarizationDecomposition()
	S:UseJonesVectorBasis()
	S:SetResolution(8)
elseif form == 'pol' then
	S:UsePolarizationDecomposition()
	S:SetResolution(8)
end

--S:SetBasisFieldDumpPrefix('basis_')
--S:SetVerbosity(1)

Pi = S:GetPoyntingFlux('Before', 0)
Pt = S:GetPoyntingFluxByOrder('After', 0)
print(S:GetNumG(), Pt[S:GetDiffractionOrder(-1,-2)][1]/Pi)

--[[
Glist = S:GetGList()
for i,v in ipairs(Glist) do
	print(i, Glist[i][1],Glist[i][2], Pt[i][1]/Pi)
end
]]

--[[
for y = -0.5*L,0.5*L,0.01*L do
	for x = -0.5*L,0.5*L,0.01*L do
		print(x,y,S:GetEpsilon({x,y,0.5}))
	end
	print('')
end
]]
