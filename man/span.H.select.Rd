% This is a template R document file
% Distributed under GPL 3 or later

\name{span.H.select}
\alias{span.H.select}
\title{A method to select the span of the sequence of lambda's for the case of heteroscedastic error}
\description{
  To compute the optimal span of the sequence of lambda's for the case of heteroscedastic error.
}
\usage{
 span.H.select(W, msigma, span=c(2,4,6,8,10,12,25), approx=FALSE, from=min(W), to=max(W), n.user=128, n.lambda=50, lambda="SJ", bw="SJ", adjust=1, na.rm = FALSE, ...)
}
\arguments{
	\item{W}{The observed data.  It is a vector of length at least 10.}
	\item{msigma}{The standard deviation \eqn{\sigma} of measurement error.  It is a single positive numeric value.}
	\item{span}{span is user-defined grids where the PDF will be evaluated.}
	\item{approx}{FALSE by default; if TRUE, the lose function uses the same one as the one of the homoscedastic error case with the average of measurement error variances}
	\item{from}{the starting point where the PDF is to be evaluated.}
	\item{to}{the starting point where the PDF is to be evaluated.}
	\item{n.user}{number of points where the PDF is to be evaluated.}
	\item{n.lambda}{number of points of lambda's.}
	\item{lambda}{Specifies the first lambda. It can be a single numeric value which has been pre-determined; or computed with the specific density bandwidth selector: 'nrd0', 'nrd', 'ucv', 'bcv', 'SJ'.}
	\item{bw}{Specifies the bandwidth of 'W'. It can be a single numeric value which has been pre-determined; or computed with the specific density bandwidth selector: 'nrd0', 'nrd', 'ucv', 'bcv', 'SJ'.}
	\item{adjust}{adjust the range there the PDF is to be evaluated.  By default, \eqn{adjust=1}.}
	\item{na.rm}{ is set to FALSE by default: no NA value is allowed.}
	\item{...}{control }
}
\details{
}
\value{
}
\author{
  X.F. Wang \email{wangx6@ccf.org}
}
\references{
Wang, X.F., Sun, J. and Fan, Z. (2011). Deconvolution density estimation with heteroscedastic errors using SIMEX.
    }
\seealso{
 \code{\link{simex.density}}.
}
\examples{
plot.simex.density <- function(X.simex,X,...){
	plot(X.simex$x, X.simex$y, type="l", xlab="x", ylab="density", lwd=3, lty=2, col=2,...)
	lines(density(X, bw="SJ"), lwd=3)
	}

############### Heteroscadestic error
N <- 1000
set.seed(123); X <- c(rnorm(N/2, mean=-2), rnorm(N/2,mean=2)); U <- rnorm(N,sd=1); msigma <- runif(N,min=0.3,max=0.5)
W <- X + msigma*U

plot.simex.density <- function(X.simex,X,...){
	plot(X.simex$x, X.simex$y, type="l", xlab="x", ylab="density", lwd=3, lty=2, col=2,...)
	lines(density(X, bw="SJ"), lwd=3)
	}


#---- Select the optimal lambda span
par(mfrow=c(1,2))
spans <- span.H.select(W, msigma, span=c(2,4,6,8,10,12,16,25), approx=TRUE)
plot(spans$span, spans$ISE, type="o", xlab="span", ylab="ISE") 

X.simex <- simex.H.density(W, msigma=msigma, adjust=1, n.lambda=50, span=spans$span[order(spans$ISE)[1]])
plot.simex.density(X.simex, X,ylim=c(0,0.25))

## more computational time needed if approx=FALSE.

# spans <- span.H.select(W, msigma, span=c(2,4,6,8,10,12,16,25), approx=FALSE)
# plot(spans$span, spans$ISE, type="o", xlab="span", ylab="ISE") 
# 
# X.simex <- simex.H.density(W, msigma=msigma, adjust=1, n.lambda=50, span=spans$span[order(spans$ISE)[1]])
# plot.simex.density(X.simex, X,ylim=c(0,0.25))


}
\keyword{smooth}

