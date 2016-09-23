-- Fig. 2b in
-- A. Christ, T. Zentgraf, J. Kuhl, S. G. Tikhodeev, N. A. Gippius, and H. Giessen
-- "Optical properties of planar metallic photonic crystal structures: Experiment and theory"
-- Physical Review B 70, 125113 (2004)

S = S4.NewSimulation()
S:SetLattice({0.5,0}, {0,0})
S:SetNumG(49)
S:AddMaterial("ITO", {2,0}) -- real and imag parts
S:AddMaterial("Quartz", {2.14,0}) -- real and imag parts
S:AddMaterial("Vacuum", {1,0})
S:AddMaterial("Gold", {-5,0})

S:AddLayer('Front', 0 , 'Vacuum')
S:AddLayer('Grating', 0.015, 'Vacuum')
S:SetLayerPatternRectangle('Grating', 'Gold', {0,0}, 0, {0.05,0})
S:AddLayer('Spacer', 0.015, 'ITO')
S:AddLayer('Back', 0, 'Quartz')
Gold = S4.GetRealMaterial('Gold')

S:SetExcitationPlanewave(
	{0,0}, -- incidence angles
	{0,0}, -- s-polarization amplitude and phase (in degrees)
	{1,0}) -- p-polarization amplitude and phase

S:UsePolarizationDecomposition()
S:SetResolution(8)

--S:OutputLayerPatternRealization('Grating', 128, 128)

Sa = {S}
npar = 3
for i = 2,npar do
	Sa[i] = S:Clone()
end

step = 0.002
for ev = 1.5,2.4,step*npar do
	for i = 1,npar do
		iev = ev+(i-1)*step
		f = 0.8065548889615557 * iev
		l2 = 1/(f*f)
		Sa[i]:SetMaterial('Gold', {Gold:GetEpsilon(iev, 'eV')})
		Sa[i]:SetMaterial('ITO', {1 + 1.81302*l2/(l2-0.07597), 0})
		Sa[i]:SetFrequency(f)
	end
	S4.SolveInParallel('Back', unpack(Sa))
	for i = 1,npar do
		iev = ev+(i-1)*step
		t = Sa[i]:GetPoyntingFlux('Back', 0)
		print(iev, -math.log(t))
	end
	io.stdout:flush()
end
