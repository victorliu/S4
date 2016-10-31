\name{S4_Layer_SetRegionHalfwidths}
\alias{S4_Layer_SetRegionHalfwidths}
\title{Set a simple region in a layer}
\description{
    Sets a region within a layer to a specified material for simple symmetric shapes.
}
\usage{
S4_Layer_SetRegionHalfwidths(S, layer, material, shape, center, halfwidths, angle_frac = 0)
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
	Can be \sQuote{interval} (for 1D only), or \sQuote{rectangle},  \sQuote{circle}, or \sQuote{ellipse} (for 2D only).
  }
  \item{center}{
    A 2-vector specifying the center of the region.
  }
  \item{halfwidths}{
    A 2-vector specifying the half-widths of the region.
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
\dontrun{S4_Layer_SetRegionHalfwidths(S, layer = Lsub, material = mSi, shape = "circle", center = c(0,0), halfwidths = c(0.2,0.2))}
}
\keyword{}
