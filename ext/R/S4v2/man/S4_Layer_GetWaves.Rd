\name{S4_Layer_GetWaves}
\alias{S4_Layer_GetWaves}
\title{Solve for planewaves}
\description{
    Solves for the set of planewaves within a specified layer.
	Planewaves may be propagating or evanescent.
}
\usage{
S4_Layer_GetWaves(S, layer)
}
\arguments{
  \item{S}{
    The Simulation object containing the layer.
  }
  \item{layer}{
    The layer in which the waves should be returned.
  }
}
\value{
    A matrix of size 11-by-(2*n), where n is the number of basis functions used.
	The 11 rows are as follows: kx, ky, kz_real, kz_imag, ux, uy, uz, cu_real, cu_imag, cv_real, cv_imag.
}
\author{
    Victor Liu
}

\examples{
\dontrun{waves = S4_Layer_GetWaves(S, Lsub)}
}
\keyword{}
