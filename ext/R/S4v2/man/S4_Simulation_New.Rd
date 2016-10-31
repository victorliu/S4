\name{S4_Simulation_New}
\alias{S4_Simulation_New}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a new Simulation object}
\description{
    Create a new Simulation object with a specified lattice and basis set.
}
\usage{
S4_Simulation_New(lattice, bases)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lattice}{
    Either a number (1D period) or a 2x2 matrix (2D lattice vectors in columns).
}
  \item{bases}{
    Either an integer (number of basis functions), a vector (basis orders for 1D), or a 2xN matrix (basis orders for 2D).
  }
}
\value{
    A Simulation object pointer.
}
\author{
    Victor Liu
}

\examples{
S <- S4_Simulation_New(lattice = matrix(c(1,0,0,1),nrow=2), bases = 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
