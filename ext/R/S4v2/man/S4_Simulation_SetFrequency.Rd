\name{S4_Simulation_SetFrequency}
\alias{S4_Simulation_SetFrequency}
\title{Set the frequency of the simulation}
\description{
    Sets the operating frequency of the simulation. When specifying wavelength, the frequency is the reciprocal of the wavelength.
}
\usage{
S4_Simulation_SetFrequency(S, freq)
}
\arguments{
  \item{S}{
    The Simulation object.
  }
  \item{freq}{
    The frequency of operation.
  }
}
\value{
    None.
}
\author{
    Victor Liu
}

\examples{
\dontrun{S4_Simulation_SetFrequency(S, 0.25)}
}
\keyword{}
