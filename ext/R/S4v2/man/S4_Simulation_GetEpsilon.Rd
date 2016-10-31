\name{S4_Simulation_GetEpsilon}
\alias{S4_Simulation_GetEpsilon}
\title{Retrieve the dielectric constant at a point in space.}
\description{
    Retrieves the dielectric constant at a specified point in space.
}
\usage{
S4_Simulation_GetEpsilon(S, p)
}
\arguments{
  \item{S}{
    The Simulation object.
  }
  \item{p}{
    The 3-vector of the position at which the dielectric constant is to be returned.
  }
}
\value{
    A complex number that is the dielectric constant at the specified point.
}
\author{
    Victor Liu
}

\examples{
\dontrun{eps = S4_Simulation_GetEpsilon(S, c(0, 0.1, 0.2))}
}
\keyword{}
