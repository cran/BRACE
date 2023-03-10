% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gendat.R
\name{gendat}
\alias{gendat}
\title{simulation data generating function}
\usage{
gendat(nsims, nobs, f, g, b, w1, w2, omegaplus, theta1, theta2, k3)
}
\arguments{
\item{nsims}{number of simulation datasets}

\item{nobs}{number of observations for one dataset}

\item{f}{parameter for generating unmeasured binary confounder}

\item{g}{parameter for generating group assignment}

\item{b}{confounder effect on group assignment}

\item{w1}{shape parameter in generating survival time for event 1
from weibull distribution}

\item{w2}{shape parameter in generating survival time for event 2
from weibull distribution}

\item{omegaplus}{multiplier on the baseline hazard for event 1}

\item{theta1}{multiplier on the baseline hazard for event 1}

\item{theta2}{multiplier on the baseline hazard for event 2}

\item{k3}{multiplier on the baseline hazard for event 2}
}
\value{
a matrix of nsims*nobs row, which consists of nsims datasets
}
\description{
generating the simulation data to apply in \code{brace}
}
\examples{
nsims = 1; nobs = 1500
f = 0.5; g = 0.333; b = 8; w1 = w2 = 0.667
theta1 = 0.5; theta2 = 1; omegaplus = 1; k3 = 0.333
sim1 = gendat(nsims,nobs,f,g,b,w1,w2,omegaplus,theta1,theta2,k3)
}
