\name{S4_Simulation_AddLayer}
\alias{S4_Simulation_AddLayer}
\title{Add a layer to a Simulation}
\description{
    Adds a new layer of a specified material and thickness to the simulation, or copy an existing layer.
	Layers may also be assigned an optional name.
	Layers are added from z = 0 towards increasing z.
}
\usage{
S4_Simulation_AddLayer(S, name, thickness, copy, material)
}
\arguments{
  \item{S}{
    The Simulation object to which the layer will be added.
  }
  \item{name}{
    An optional name for the layer.
  }
  \item{thickness}{
    The thickness of the layer.
  }
  \item{copy}{
    A Layer object to copy. When specifying \code{copy}, \code{material} should not be specified.
  }
  \item{material}{
    A Material object specifying the base material of the layer. When specifying \code{material}, \code{copy} should not be specified.
  }
}
\value{
    A Layer object pointer.
}
\author{
    Victor Liu
}

\examples{
\dontrun{Lsub <- S4_Simulation_AddLayer(S, name = "Substrate", thickness = 0.7, material = mSi)}
}
\keyword{}
