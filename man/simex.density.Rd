% This is a template R document file
% Distributed under GPL 3 or later

\name{simex.density}
\alias{simex.density}
\title{Estimating probability density function from data with homoscedastic measurement error using SIMEX}
\description{
  To compute the probability density function from data contaminated with homoscedastic measurement error using SIMEX.
}
\usage{
simex.density(W, msigma, x, from, to, n.user=128, n.lambda=50, 
	lambda="SJ", span=8, adjust=1, na.rm = FALSE, ...)
}
\arguments{
  \item{W}{The observed data.  It is a vector of length at least 10.}
  \item{msigma}{The standard deviation \eqn{\sigma} of measurement error.  It is a single positive numeric value.}
  \item{x}{x is user-defined grids where the PDF will be evaluated.}
  \item{from}{the starting point where the PDF is to be evaluated.}
  \item{to}{the starting point where the PDF is to be evaluated.}
  \item{n.user}{number of points where the PDF is to be evaluated.}
  \item{n.lambda}{number of points of lambda's.}
  \item{lambda}{Specifies the first lambda. It can be a single numeric value which has been pre-determined; or computed with the specific density bandwidth selector: 'nrd0', 'nrd', 'ucv', 'bcv', 'SJ'.}
  \item{span}{Specifies the span of lambda.}
  \item{adjust}{adjust the range there the PDF is to be evaluated.  By default, \eqn{adjust=1}.}
  \item{na.rm}{ is set to FALSE by default: no NA value is allowed.}
  \item{...}{control }
}
\details{
 }
\value{
An object of class ``NPsimex''.
}
\author{
  X.F. Wang \email{wangx6@ccf.org}
}
\references{
Wang, X.F., Sun, J. and Fan, Z. (2011). Deconvolution density estimation with heteroscedastic errors using SIMEX.
}
\seealso{
 \code{\link{simex.density}}, \code{\link{lambda.select}}, \code{\link{span.select}}.
}
\examples{

############### Homoscadestic error
N <- 1000
set.seed(123); X <- c(rnorm(N/2, mean=-2), rnorm(N/2,mean=2)); U <- rnorm(N,sd=1)
msigma <- 0.5
W <- X + msigma*U

plot.simex.density <- function(X.simex,X,...){
	plot(X.simex$x, X.simex$y, type="l", xlab="x", ylab="density", lwd=3, lty=2, col=2,...)
	lines(density(X, bw="SJ"), lwd=3)
	}

par(mfrow=c(2,2))
X.simex1 <- simex.density(W, msigma=msigma, adjust=1, n.lambda=50, span=1)
plot.simex.density(X.simex1, X, ylim=c(0,0.25))

X.simex2 <- simex.density(W, msigma=msigma, adjust=1, n.lambda=50, span=3)
plot.simex.density(X.simex2, X, ylim=c(0,0.25))

X.simex3 <- simex.density(W, msigma=msigma, adjust=1, n.lambda=50, span=8)
plot.simex.density(X.simex3, X, ylim=c(0,0.25))

X.simex4 <- simex.density(W, msigma=msigma, adjust=1, n.lambda=50, span=35)
plot.simex.density(X.simex4, X, ylim=c(0,0.25))

}
\keyword{nonparametric}
\keyword{smooth}
\keyword{density estimation}
\keyword{measurement error}

