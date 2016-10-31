\name{S4_Layer_SetRegionVertices}
\alias{S4_Layer_SetRegionVertices}
\title{Set a region in a layer defined by vertices}
\description{
    Sets a region within a layer to a specified material for a shape defined by vertices.
}
\usage{
S4_Layer_SetRegionVertices(S, layer, material, shape, center, vertices, angle_frac = 0)
}
\arguments{
  \item{S}{
    The Simulation object that contains the layer.
  }
  \item{layer}{
    The Layer object to modify.
  }
  \item{material}{
    The Material object to set the region to.
  }
  \item{shape}{
    A string specifying the shape of the region.
	Can be \sQuote{polygon}.
  }
  \item{center}{
    A 2-vector specifying the center of the region.
  }
  \item{vertices}{
    A 2-by-n matrix specifying the vertices of the shape, where each vertex is one column.
	The vertices should be oriented in counter-clockwise order around the shape.
  }
  \item{angle_frac}{
    A number specifying the rotation of the region, from 0 to 1 for one rotation.
  }
}
\value{
    None.
}
\author{
    Victor Liu
}

\examples{
\dontrun{S4_Layer_SetRegionVertices(S, layer = Lsub, material = mSi, shape = "circle", center = c(0,0), vertices = matrix(c(0,0, 1,0, 0,1), nrow=2))}
}
\keyword{}
