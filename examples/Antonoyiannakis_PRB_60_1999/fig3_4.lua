-- Reproduces Fig. 3 and 4 in
-- M. I. Antonoyiannakis and J. B. Pendry,
-- "Electromagnetic forces in photonic crystals",
-- Phys. Rev. B, Vol. 60, No. 4, 1999.

S = S4.NewSimulation()
S:SetLattice({1,0}, {0,1})
S:SetNumG(1)
S:AddMaterial('Al', {1,0}) -- real and imag parts
S:AddMaterial('Vacuum', {1,0})

S:AddLayer('Above', 0 , 'Vacuum')
S:AddLayer('Below', 0, 'Al')

S:SetExcitationPlanewave(
	{0,0}, -- incidence angles
	{1,0}, -- s-polarization amplitude and phase (in degrees)
	{0,0}) -- p-polarization amplitude and phase

f_p = 15
gamma = 0.1

for f=1,25,0.2 do
	S:SetFrequency(f)
	-- eps = 1 - f_p^2/(f*(f + i*gamma))
	-- eps = 1 - f_p^2(f - i*gamma)/(f*(f^2 + gamma^2))
	epsr = 1 - f_p*f_p/(f*f + gamma*gamma)
	epsi = f_p*f_p/(f*f + gamma*gamma) * gamma/f
	S:SetMaterial('Al', {epsr,epsi});
	
	-- Reflected power
	forward,reflected = S:GetPoyntingFlux('Above', 0)
	reflected = -reflected/forward
	
	T1x,T1y,T1z = S:GetStressTensorIntegral('Above', 0)
	T2x,T2y,T2z = S:GetStressTensorIntegral('Below', 0)
	
	-- We need to plot -T1z because the surface normal is in the -z direction.
	print (f .. '\t' .. reflected .. '\t' .. -T1z .. '\t' .. T2z)
end
