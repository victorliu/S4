\name{S4_Simulation_AddMaterial}
\alias{S4_Simulation_AddMaterial}
\title{Add a material to a Simulation}
\description{
    Adds a new material with a specified dielectric constant, and optionally, a name.
	The dielectric constant may be a scalar real or complex value, or a 3x3 real or complex tensor.
}
\usage{
S4_Simulation_AddMaterial(S, name, epsilon)
}
\arguments{
  \item{S}{
    The Simulation object to which the material will be added.
  }
  \item{name}{
    An optional name for the material.
  }
  \item{epsilon}{
    The dielectric constant of the material. This value may be a scalar or a 3x3 matrix.
  }
}
\value{
    A Material object pointer.
}
\author{
    Victor Liu
}

\examples{
\dontrun{mSi <- S4_Simulation_AddMaterial(S, name = "Silicon", epsilon = 12)}
}
\keyword{}
