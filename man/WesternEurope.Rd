\name{WesternEurope}
\alias{WesternEurope_marine}
\alias{WesternEurope_terrestrial}
\alias{WesternEurope}
\docType{package}
\title{
Western European rock outcrop bias}
\description{
Quantification of the Western European rock outcrop of the Phanerozoic
}
\usage{
data(WesternEurope_marine)
data(WesternEurope_terrestrial)
data(WesternEurope)
}
\format{
A data frame with 72 values for preservation rate bias (relr) and the younger end of the time interval (py) for which this preservation rate has been calculated.
}
\details{
Western European marine and terrestrial rock outcrops of the Phanerozoic have been quantified by Smith & McGowan (2007) and McGowan & Smith (2008) as the number of geological maps of France, Spain, England, and Wales containing this type of rock outcrop of a given time interval (using bins with an average size of 7.5 myr). This quantification is here translated into maps per myr and is scaled so that the maximum number of maps per myr equals 1. This data serves as a proxy for preservation rate bias. It is used by functions adjust.pres.rate and find.prior when called with corr="WesternEurope_marine", corr="WesternEurope_terrestrial", or corr="WesternEurope" (which combines both the marine and terrestrial rock outcrop).
}
\author{
Michael Matschiner

Maintainer: Michael Matschiner <michaelmatschiner@mac.com>
}
\source{
Matschiner M (2011) in prep.
}
\references{
McGowan AJ, Smith AB (2008) \emph{Paleobiology} 34:80-103.\cr
Smith AB, McGowan AJ (2007) \emph{Palaeontology} 50:765-774.
}
\keyword{datasets}
