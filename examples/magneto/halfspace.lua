-- Infinite halfspace of magneto-optic material.
-- Assume the permittivity tensor in the receiving halfspace is
--   [  e   ie' 0 ]
--   [ -ie' e   0 ]
--   [  0   0   e ]
-- Then solving the wave equation gives two circular polarization modes:
--   E = (1, i)/sqrt(2), km = omega sqrt(e-ie)
--   E = (1,-i)/sqrt(2), kp = omega sqrt(e+ie)
-- For normal incidence, we match boundary conditions:
--   (E0+Erx, Ery) = a(1,i)/sqrt(2) + b(1,-i)/sqrt(2)
--   k0(E0-Erx, -Ery) = a km(1,i)/sqrt(2) + b kp(1,-i)/sqrt(2)
--   k0 = omega
-- Solving these gives
--   a = E0 k0/(k0+km)/sqrt(2), b = E0 k0/(k0+kp)/sqrt(2)
--   Er = E0/((k0+km)(k0+kp))/sqrt(2) (k0 k0 - km kp, i k0(kp-km))/sqrt(2)
-- The transmitted field is
--   Etx = Re[   E0 k0/sqrt(2) (exp(i km z)/(k0+km) + exp(i kp z)/(k0+kp)) ]
--   Ety = Re[ i E0 k0/sqrt(2) (exp(i km z)/(k0+km) - exp(i kp z)/(k0+kp)) ]


S = S4.NewSimulation()
S:SetLattice({1,0}, {0,1})
S:SetNumG(1)
f = 0.4
S:SetFrequency(f)

eps  = 10
epsp = 3
S:AddMaterial("Dielectric", {
	{eps, 0}, {0, epsp}, {0, 0},
	{0, -epsp}, {eps, 0}, {0, 0},
	{0, 0}, {0, 0}, {eps, 0}
	})
S:AddMaterial("Vacuum", {1,0})

S:AddLayer('AirAbove', 0, 'Vacuum')
S:AddLayer('Stuff', 0, 'Dielectric')

S:SetExcitationPlanewave(
	{0,0}, -- incidence angles (spherical coordinates: phi in [0,180], theta in [0,360])
	{0,0}, -- s-polarization amplitude and phase (in degrees)
	{1,0}) -- p-polarization amplitude and phase

k0 = 2*math.pi*f
km = k0*math.sqrt(eps-epsp)
kp = k0*math.sqrt(eps+epsp)

for z = 0,10,0.01 do
	Exr,Eyr = S:GetEField({0,0,z})
	Etx =  (math.cos(km*z)/(k0+km) + math.cos(kp*z)/(k0+kp)) / f
	Ety = -(math.sin(km*z)/(k0+km) - math.sin(kp*z)/(k0+kp)) / f
	print(z, Exr,Etx, Eyr,Ety)
end
