\name{adjust.pres.rate}
\alias{adjust.pres.rate}
\docType{package}
\title{
Preservation rates with rock outcrop bias}
\description{
Adjusts the preservation rate over a given period to the rock outcrop bias
}
\usage{
adjust.pres.rate(r,ty,t1,t2,corr,plot)
}
\arguments{
\item{r}{
Preservation rate.
}
\item{ty}{
The younger age boundary of the geological formation in which the oldest fossil constraining the node has been found (in Ma). This will be used as the prior offset.
}
\item{t1}{
Younger end of the range for which the preservation rate has been calculated.
}
\item{t2}{
Older end of the range for which the preservation rate has been calculated.
}
\item{corr}{
The applied correction for outcrop bias. Options include:\cr
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
"Australia_terrestrial" (McGowan & Smith 2008), and\cr
"Australia" (McGowan & Smith 2008).\cr
When using the adjusted preservation rate in function find.prior, the same correction must be applied there.
}
\item{plot}{
This optional parameter allows to suppress graphical output with plot=FALSE. The default setting is plot=TRUE.
}
}
\details{
The function find.prior allows for adjustment of prior distributions to bias in the fossil record that results from the
heterogeneous availability of sedimentary rock outcrops per interval. If a correction is applied in find.prior, the 
preservation rate should also be adjusted to any rock outcrop bias at the time of the clade-constraining fossil. For
example, if a preservation rate of 0.03 has been calculated for Cenozoic mammals, then this specifies the average rate
between 65.5-0 Ma. If the oldest fossil of a given mammal clade dates to the Eocene (56-34 Ma), then find.prior would
calculate an lognormal, gamma, or exponential prior distribution with offset ty=34 Ma, based on a parameter combination
including the parameters r (preservation rate) and p_q (p-q, the net diversification rate). If however, the preservation
rate at ty=34 Ma can be expected to be higher or lower than the average preservation rate (because higher- or lower-than-
average amounts of rock outcrops are known of this age), then the preservation rate used to define the resulting prior
distribution should be adjusted. Say, there is twice as much rock outcrop of an age of ty=34 Ma compared to the average
outcrop amount for the whole Cenozoic, then the preservation rate used to define any prior with offset ty=34 should be
twice the overall Cenozoic preservation rate. In the example given above, it should be set to 0.06.
The function adjust.pres.rate facilitates this calculation. Given parameters r (preservation rate - 0.03 in the above
example), ty (prior offset / youngest possible age of the clade's oldest fossil, 34 Ma in the example), t1 (younger end
of the range for which the preservation rate has been calculated - 0 Ma in the example of Cenozoic mammals, as the Cenozoic
extends to the present), t2 (older end of the range for which the preservation rate has been calculated - 65.5 Ma, the
begin of the Cenozoic in the example), and corr (the applied correction for outcrop bias, in the example, it could be 
"global_terrestrial"), the function adjust.pres.rate calculates a corrected preservation rate that should then be used
as a parameter in the function \code{\link{find.prior}}.
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
Smith AB, McGowan AJ (2007) \emph{Palaeontology} 50:765-774.\cr
Wall PD, Ivany LC, Wilkinson BH (2009) \emph{Paleobiology} 35:146-167.
}
\keyword{methods}
\seealso{
\code{\link{find.prior}}
}
\examples{
adjust.pres.rate(r=0.02,ty=34,t1=0,t2=65.5,corr="global_terrestrial")
}
