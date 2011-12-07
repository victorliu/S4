-- Fig. 4 in
-- Quasi-guided modes and optical properties of photonic crystal slabs
-- S.G. Tikhodeev, A.L. Yablonskii, E.A. Muljarov, N.A. Gippius, and T. Ishihara
-- Phys. Rev. B 66, 45102 (2002)

S = S4.NewSimulation()
S:SetLattice({0.68,0}, {0,0.68})
S:SetNumG(100)
S:AddMaterial("Active", {3.97,0}) -- real and imag parts
S:AddMaterial("Quartz", {2.132,0}) -- real and imag parts
S:AddMaterial("Vacuum", {1,0})

S:AddLayer('Front', 0 , 'Vacuum')
S:AddLayer('Slab', 0.12, 'Quartz')
S:SetLayerPatternRectangle('Slab', 'Active', {0,0}, 0, {0.4*0.68,0.4*0.68})
S:AddLayer('Back', 0, 'Quartz')

S:SetExcitationPlanewave(
	{0,0}, -- incidence angles
	{1,0}, -- s-polarization amplitude and phase (in degrees)
	{0,0}) -- p-polarization amplitude and phase

S:UsePolarizationDecomposition()
S:SetResolution(8)

--S:OutputLayerPatternRealization('Slab', 128, 128)

Sa = {S}
npar = 3
for i = 2,npar do
	Sa[i] = S:Clone()
end

step = 0.002
for ev = 1,2.6,step*npar do
	for i = 1,npar do
		f = 0.8065548889615557 *(ev+(i-1)*step)
		Sa[i]:SetFrequency(f)
	end
	S4.SolveInParallel('Back', unpack(Sa))
	for i = 1,npar do
		t = Sa[i]:GetPoyntingFlux('Back', 0)
		print(ev+(i-1)*step, t)
	end
	io.stdout:flush()
end
