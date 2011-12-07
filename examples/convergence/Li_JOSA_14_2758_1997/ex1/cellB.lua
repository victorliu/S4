pcall(loadstring(S4.arg))

if not NG then NG = 100 end

sidelen = 1.25
L = sidelen * math.sqrt(2)

local S = S4.NewSimulation()
S:SetLattice({L,0}, {0,L})
S:SetNumG(NG)
S:AddMaterial("Dielectric", {1.5*1.5,0}) -- real and imag parts
S:AddMaterial("Vacuum", {1,0})

S:AddLayer('Before', 0 , 'Dielectric')
S:AddLayer('Slab', 1, 'Vacuum')
S:SetLayerPatternRectangle('Slab', 'Dielectric', {0,0}, 45, {0.5*sidelen, 0.5*sidelen})
S:AddLayer('After', 0 , 'Vacuum')
S:SetFrequency(1)

S:SetExcitationPlanewave(
	{0,45}, -- incidence angles
	{1,0}, -- s-polarization amplitude and phase (in degrees)
	{0,0}) -- p-polarization amplitude and phase

if form == 'ref' then
	S:UseDiscretizedEpsilon(false)
elseif form == 'fft' then
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
print(S:GetNumG(), Pt[2][1]/Pi)


--[[
for y = -0.5*L,0.5*L,0.01*L do
	for x = -0.5*L,0.5*L,0.01*L do
		print(x,y,S:GetEpsilon({x,y,0.5}))
	end
	print('')
end
]]
