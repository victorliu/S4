\name{S4v2-package}
\alias{S4v2-package}
\docType{package}
\title{
R interface to S4
}
\description{
This package provides an R interface to S4, a package for solving Maxwell's equations using the method of Rigorous Coupled Wave Analysis.
}
\details{
\tabular{ll}{
Package: \tab S4v2\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-10-29\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Victor Liu

Maintainer: Victor Liu <victor.liu@gmail.com>
}
\keyword{ package }
\examples{
S <- S4_Simulation_New(lattice = matrix(c(1,0,0,1),nrow=2), bases = 100)

mVac = S4_Simulation_AddMaterial(S, name = "Vacuum", epsilon = 1);
mSi  = S4_Simulation_AddMaterial(S, name = "Silicon", epsilon = 12);

L0 = S4_Simulation_AddLayer(S, name = "AirAbove", thickness = 0, material = mVac);
Ls = S4_Simulation_AddLayer(S, name = "Slab", thickness = 0.5, material = mSi);
S4_Layer_SetRegionHalfwidths(S, layer = Ls, shape = "circle", material = mVac, halfwidths = c(0.2, 0.2));
LN = S4_Simulation_AddLayer(S, name = "AirBelow", thickness = 0, copy = L0);

S4_Simulation_ExcitationPlanewave(S, k = c(0,0,1), u = c(1,0,0), cu = 1, cv = 0);

\dontrun{
for(ifreq in 250:600){
    freq = ifreq/1000
    S4_Simulation_SetFrequency(S, freq);
    p = S4_Layer_GetPowerFlux(S, layer = L0);
    power_inc = Re(p[1]);
    power_ref = Re(p[2]);
    p = S4_Layer_GetPowerFlux(S, layer = LN);
    power_tra = Re(p[1]);
    cat(sprintf("%f\t%f\t%f\n", freq, power_tra/power_inc, power_ref/power_inc))
}
}

}
