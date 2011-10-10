\name{NorthAmerica}
\alias{NorthAmerica_marine}
\alias{NorthAmerica_terrestrial}
\alias{NorthAmerica}
\docType{package}
\title{
North American rock outcrop bias}
\description{
Quantification of the North American rock outcrop of the Phanerozoic
}
\usage{
data(NorthAmerica_marine)
data(NorthAmerica_terrestrial)
data(NorthAmerica)
}
\format{
Data frames with 23 values for preservation rate bias (relr) and the younger end of the time interval (py) for which this preservation rate has been calculated.
}
\details{
North American marine and terrestrial rock outcrops of the Phanerozoic have been quantified by Wall et al. (2009) by combination of data on North American outcrop areas from the UNESCO \emph{Geological Atlas of the World} (Choubert & Faure-Muret 1976) with volumetric data on North American sediment compositions (Ronov 1980, 1994). Bin size for the 22 used phanerozoic time intervals varies between 6 (Early Triassic) and 54 (Cambrian) myr, with an average length of 24.4 myr. The 23rd value results from an extrapolation of Miocene rock outcrop to the present. Rock outcrop quantification is here translated into area per myr and is scaled so that the maximum area per myr equals 1. This data serves as a proxy for preservation rate bias. It is used by functions adjust.pres.rate and find.prior when called with corr="NorthAmerica_marine", corr="NorthAmerica_terrestrial", or corr="NorthAmerica" (which combines both the marine and terrestrial rock outcrop).
}
\author{
Michael Matschiner

Maintainer: Michael Matschiner <michaelmatschiner@mac.com>
}
\source{
Matschiner M (2011) in prep.
}
\references{
Choubert G, Faure-Muret A (1976) \emph{Geological atlas of the world. 1:10,000,000. 22 sheets with explanations}. UNESCO Commission for the Geological Map of the World, Paris.\cr
Ronov AB (1980) The earth's sedimentary shell: quantitative patterns of its structures, compositions and evolution (the 20th V. I. Vernadskiy Lecture). pp. 1-80 in Yaroshevskiy AA (ed.) \emph{The earth's sedimentary shell}. Nauka, Moscow. American Geological Institute Reprint Series 5:1-73.\cr
Ronov AB (1994) Phanerozoic transgressions and regressions on the continents: a quantitative approach based on areas flooded by the sea and areas of marine and continental deposition. \emph{American Journal of Science} 294:777-801.\cr
Wall PD, Ivany LC, Wilkinson BH (2009) \emph{Paleobiology} 35:146-167.
}
\keyword{datasets}
