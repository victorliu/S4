-- One curve of Fig. 2b in
-- Wonjoo Suh, M. F. Yanik, Olav Solgaard, and Shanhui Fan,
-- "Displacement-sensitive photonic crystal structures based on guided resonance in photonic crystal slabs",
-- Appl. Phys. Letters, Vol. 82, No. 13, 2003

S = S4.NewSimulation()
S:SetLattice({1,0}, {0,1})
S:SetNumG(40)
S:AddMaterial("Silicon", {12,0}) -- real and imag parts
S:AddMaterial("Vacuum", {1,0})

S:AddLayer('AirAbove', 0 , 'Vacuum')
S:AddLayer('Slab', 0.55, 'Silicon')
S:SetLayerPatternCircle('Slab', 'Vacuum', {0,0}, 0.4)
S:AddLayerCopy('Spacer', 1.1, 'AirAbove')
S:AddLayerCopy('Slab2', 0.55, 'Slab')
S:AddLayerCopy('AirBelow', 0, 'AirAbove')

S:SetExcitationPlanewave(
	{0,0}, -- incidence angles
	{1,0}, -- s-polarization amplitude and phase (in degrees)
	{0,0}) -- p-polarization amplitude and phase

S:UsePolarizationDecomposition()

for freq=0.49,0.6,0.001 do
		S:SetFrequency(freq)
		forward = S:GetPoyntingFlux('AirBelow', 0)
		print (freq .. '\t' .. forward)
end
