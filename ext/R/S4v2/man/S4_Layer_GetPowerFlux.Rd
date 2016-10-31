\name{S4_Layer_GetPowerFlux}
\alias{S4_Layer_GetPowerFlux}
\title{Solve for forward and backward fluxes}
\description{
    Solves for the total forward and backward (in terms of z) power flux in a specified layer.
}
\usage{
S4_Layer_GetPowerFlux(S, layer, offset = 0)
}
\arguments{
  \item{S}{
    The Simulation object containing the layer.
  }
  \item{layer}{
    The layer in which the fluxes should be returned.
  }
  \item{offset}{
    An optional offset into the layer at which the fluxes should be computed.
	The default value is 0.
  }
}
\value{
    A 2-vector of the complex forward and backward propagating powers in the layer.
}
\author{
    Victor Liu
}

\examples{
\dontrun{powers = S4_Layer_GetPowerFlux(S, Lsub)}
}
\keyword{}
