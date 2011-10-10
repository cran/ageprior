\name{net.div.rate}
\alias{net.div.rate}
\docType{package}
\title{
Exponential and saturated net diversification rate}
\description{
Calculates net diversification rates assuming a model of exponential diversification, with or without saturation
}
\usage{
net.div.rate(n,dt,sat)
}
\arguments{
\item{n}{
Extant species diversity of the investigated clade.
}
\item{dt}{
Time interval over which the extant species diversity is believed to have accumulated (= the clade age).}
\item{sat}{
Saturation (optional). If set to TRUE, then it is assumed that half the current species richness existed after half the clade's age, and that diversity slowly saturated after this. This is similar to logistic growth. By default, sat is set to FALSE.
}
}
\details{
The function net.div.rate calculated the net diversification rate for a clade, given n (the number of extant species in this clade), dt (the time since clade origin), and assuming an exponential diversification model. Optionally, it can be assumed that half the extant species richness was present after half the time, and that since that time, diversity has saturated. In this case, the parameter sat should be set to TRUE. Estimation of the net diversification is problematic, because it requires an estimate for the age of the investigated clade, which one typically tries to find with this type of analysis. Thus, the age estimate influences the net diversification rate estimate, which in turn again influences the age estimate. I recommend to estimate the net diversification rate before and after the molecular dating analysis, and repeat the analysis if the a posteriori net diversification rate is too different from the a priori assumptions. It may be very useful to run an analysis with the software MEDUSA (Alfaro et al. 2009), or the forthcoming FOSSILMEDUSA (www.webpages.uidaho.edu/~lukeh/software/index.html), to gain insights into the net diversification rate dynamics of the investigated clade. Experience has shown that molecular dating results are relatively robust to the chosen net diversification rate (higher estimates of net diversification lead to older age estimates). It may be wise to try a range of possible net diversification rates.
}
\author{
Michael Matschiner

Maintainer: Michael Matschiner <michaelmatschiner@mac.com>
}
\source{
Matschiner M (2011) in prep.
}
\references{
Alfaro ME, Santini F, Brock CD, et al. (2009) \emph{Proc Natl Acad Sci} 106:13410-13414.
}
\keyword{methods}
\seealso{
\code{\link{find.prior}}
}
\examples{
net.div.rate(n=4500,dt=180,sat=TRUE)
}
