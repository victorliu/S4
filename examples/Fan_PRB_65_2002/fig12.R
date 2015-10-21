# Bottom pane of Fig. 12 in
# Shanhui Fan and J. D. Joannopoulos,
# "Analysis of guided resonances in photonic crystal slabs",
# Phys. Rev. B, Vol. 65, 235112

require("RCWA")
S <- new(Simulation,matrix(c(1,0,0,1),nrow=2,ncol=2), 100)

Mvac = S$set_material("Vacuum", 1)
MSi  = S$set_material("Silicon", 12 + 0i)

Ltop  = S$add_layer("AirAbove", 0 , Mvac)
Lslab = S$add_layer("Slab", 0.5, MSi)
        S$pattern_circle(Lslab, Mvac, c(0,0), 0.2)
Lbot  = S$add_layer("AirBelow", 0, Mvac)

S$excitation_planewave(
	c(0,0), # incidence angles (radians)
	1, # s-polarization (complex amplitude)
	0) # p-polarization (complex amplitude)

for (freq in seq(0.25, 0.6, length.out=100)){
	S$set_frequency(freq)
	powers_top = S$get_power_flux(Ltop, 0)
	powers_bot = S$get_power_flux(Lbot, 0)
	cat(sprintf("%f\t%f\t%f\n", freq, powers_bot[["forward"]], powers_top[["backward"]]))
}
