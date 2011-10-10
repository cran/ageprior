\name{Australia}
\alias{Australia_marine}
\alias{Australia_terrestrial}
\alias{Australia}
\docType{package}
\title{
Australian rock outcrop bias}
\description{
Quantification of the Australian rock outcrop of the Phanerozoic
}
\usage{
data(Australia_marine)
data(Australia_terrestrial)
data(Australia)
}
\format{
Data frames with 72 values for preservation rate bias (relr) and the younger end of the time interval (py) for which this preservation rate has been calculated.
}
\details{
Australian marine and terrestrial rock outcrops of the Phanerozoic have been quantified by McGowan & Smith (2008) as the number of geological maps of Australia containing this type of rock outcrop of a given time interval (using bins with an average size of 7.5 myr). This quantification is here translated into maps per myr and is scaled so that the maximum number of maps per myr equals 1. This data serves as a proxy for preservation rate bias. It is used by functions adjust.pres.rate and find.prior when called with corr="Australia_marine", corr="Australia_terrestrial", or corr="Australia" (which combines both the marine and terrestrial rock outcrop).
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
}
\keyword{datasets}
