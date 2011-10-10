\name{find.prior}
\alias{find.prior}
\docType{package}
\title{
Prior distributions for molecular dating}
\description{
Calculation of the parameters of Bayesian prior distributions for divergence dates based on fossil calibrations
}
\usage{
find.prior(id,ty,to,r,p_q,lag,corr,plot)
}
\arguments{
\item{id}{
A name for the node. If this matches the name of a taxonomic group defined in BEAUTi, then the output of this script can be integrated directly into the BEAST XML file to constrain the tmrca of this group.
}
\item{ty}{
The younger age boundary of the geological formation in which the oldest fossil constraining the node has been found (in Ma). This will be used as the prior offset.
}
\item{to}{
The older age boundary of the geological formation in which the oldest fossil constraining the node has been found (in Ma).
}
\item{r}{
Preservation rate (see Foote 1997, Foote et al. 1999). If a correction for outcrop bias is applied (see below), the preservation rate may need to be adjusted with \code{\link{adjust.pres.rate}} before running find.prior.
}
\item{p_q}{
An estimate for the net diversification rate p-q (origination rate - extinction rate) of the investigated clade. This can be calculated from the number of extant species, an age estimate, and with assumption of an exponential diversification model (with or without saturation) using the function \code{\link{net.div.rate}}.
}
\item{lag}{
An optional lag time parameter to account for the facts that directly after speciation fossil preservation is less likely due to potentially smaller population sizes and geographic ranges, and for the decreased probability of correct taxonomic assignment of fossils that preserved shortly after the origin of a given clade. The lag parameter specifies the mean of an exponential distribution for the delay between clade origin and the point in time, from which on the above effects become negligible with respect to the preservation probability. These effects can be ignored completely by specifying lag=0 (the default).
}
\item{corr}{
Allows an optional correction for bias in the fossil record that results from the heterogeneous availability of sedimentary rock outcrops per interval. Options include:\cr
"global_marine" (Wall et al. 2009),\cr
"global_terrestrial" (Wall et al. 2009),\cr
"global" (Wall et al. 2009),\cr
"NorthAmerica_marine" (Wall et al. 2009),\cr
"NorthAmerica_terrestrial" (Wall et al. 2009),\cr
"NorthAmerica" (Wall et al. 2009),\cr
"WesternEurope_marine" (Smith & McGowan 2007),\cr
"WesternEurope_terrestrial" (Smith & McGowan 2007),\cr
"WesternEurope" (Smith & McGowan 2007),\cr
"Australia_marine" (McGowan & Smith 2008),\cr
"Australia_terrestrial" (McGowan & Smith 2008),\cr
"Australia" (McGowan & Smith 2008), and\cr
"none".\cr
Unless this is set to "none", care should be taken that the preservation rate is adjusted to rock outcrop bias at the youngest possible age of the oldest fossil of a given clade. This adjustment is performed by function \code{\link{adjust.pres.rate}}.
}
\item{plot}{
This optional parameter allows to suppress graphical output with plot=FALSE. The default setting is to plot the calculated probability distribution for linage nonpreservation, and the best log-normal, gamma, or exponential distribution approximation.}
}
\details{
This function finds the best log-normal, gamma, or exponential approximation to probabilistic bounds on lineage nonpreservation according to Foote et al. (1999). Offset and mean of the log-normal, gamma, or exponential distribution are fixed to match those of the theoretical distribution Rt (Matschiner 2011), and so is the offset. As a measure of fit between distributions, the root mean square deviation is calculated to decide which type of distribution provides the best approximation. For log-normal distributions parameter sigma is chosen to minimize the root mean square deviation, and the same is performed for the shape parameter of gamma distributions. The three types of distributions are readily implemented into divergence date estimation with the software BEAST (Drummond & Rambaut 2007), and output is written that can be inserted into BEAST XML files to apply the chosen prior distributions.
}
\author{
Michael Matschiner

Maintainer: Michael Matschiner <michaelmatschiner@mac.com>
}
\source{
Matschiner M (2011) in prep.
}
\references{
Drummond AJ, Rambaut A (2007) \emph{BMC Evol Biol} 7:214.\cr
Foote M (1997) \emph{Paleobiology} 23:278-300.\cr
Foote M, Hunter JP, Janis CM, Sepkoski Jr J (1999) \emph{Science} 283:1310-1314.\cr
McGowan AJ, Smith AB (2008) \emph{Paleobiology} 34:80-103.\cr
Smith AB, McGowan AJ (2007) \emph{Palaeontology} 50:765-774.\cr
Wall PD, Ivany LC, Wilkinson BH (2009) \emph{Paleobiology} 35:146-167.
}
\keyword{methods}
\seealso{
\code{\link{adjust.pres.rate}}
\code{\link{net.div.rate}}
}
\examples{
find.prior(id="NodeA",ty=10,to=20,r=0.02,p_q=0.02,corr="global_marine")
find.prior(id="NodeB",ty=65.5,to=70.6,r=0.02,p_q=0.01,lag=2)
find.prior(id="NodeC",ty=145.5,to=150.8,r=0.02,p_q=0.04)
}
