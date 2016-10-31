\name{S4_Simulation_ExcitationPlanewave}
\alias{S4_Simulation_ExcitationPlanewave}
\title{Set a planewave excitation condition}
\description{
    Sets the parameters for an incident planewave condition.
	A planewave is specified by its k-vector, a polarization vector (u, orthogonal to k), and complex amplitudes of the electric field along u and v (where u cross v is in the direction of k).
}
\usage{
S4_Simulation_ExcitationPlanewave(S, k, u, cu, cv)
}
\arguments{
  \item{S}{
    The Simulation object.
  }
  \item{k}{
    The 3-vector of the direction of the planewave's k-vector.
	This vector does not need to be normalized.
  }
  \item{u}{
    The 3-vector of the direction of the planewave's u polarization vector.
	This vector does not need to be normalized, and will be orthogonalized against k.
	The v direction is implicity defined as that for which u cross v is in the direction of k.
  }
  \item{cu}{
	The complex amplitude of the electric field along the u direction.
  }
  \item{cv}{
	The complex amplitude of the electric field along the v direction.
  }
}
\value{
    None.
}
\author{
    Victor Liu
}

\examples{
\dontrun{S4_Simulation_ExcitationPlanewave(S, k = c(0,0,1), u = c(1,0,0), cu = 1, cv = 1i)}
}
\keyword{}
