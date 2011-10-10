\name{ageprior}
\alias{ageprior}
\docType{package}
\title{
Prior distributions for molecular dating
}
\description{
Calculation of the parameters of Bayesian prior distributions for divergence dates based on fossil calibrations. Necessary parameters are preservation rate, net diversification rate, and the age of the oldest fossil of the investigated clade. The preservation rate can be calculated from the fossil record (see Foote 1997, Foote et al. 1999), and should be adjusted with \code{\link{adjust.pres.rate}} if a correction for rock outcrop bias (option "corr" of function \code{\link{find.prior}}) is used. The net diversification rate can be estimated from the number of extant taxa and a rough age estimate for the clade (assuming an exponential diversification model with or without saturation) using \code{\link{net.div.rate}}, or with more sophisticated models implemented in MEDUSA (Alfaro et al. 2009; part of the R package Geiger). Once these parameter estimates are at hand, they can be used to calculate an age prior distribution with \code{\link{find.prior}}.
}
\details{
\tabular{ll}{
Package: \tab ageprior\cr
Type: \tab Package\cr
Version: \tab 0.92\cr
Date: \tab 2011-10-10\cr
License: \tab GPL >= 2\cr
LazyLoad: \tab yes\cr
 }
}

\author{
Michael Matschiner

Maintainer: Michael Matschiner <michaelmatschiner@mac.com>
}

\source{
Matschiner M (2011) in prep.
}
\references{
Alfaro ME, Santini F, Brock CD, et al. (2009) \emph{Proc Natl Acad Sci} 106:13410-13414.\cr
Foote M (1997) \emph{Paleobiology} 23:278-300.\cr
Foote M, Hunter JP, Janis CM, Sepkoski Jr J (1999) \emph{Science} 283:1310-1314.
}
\keyword{package}
\seealso{
\code{\link{find.prior}}
\code{\link{adjust.pres.rate}}
\code{\link{net.div.rate}}
}
