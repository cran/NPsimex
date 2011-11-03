% This is a template R document file
% Distributed under GPL 3 or later

\name{lambda.select}
\alias{lambda.select}
\title{A rule of thumb lambda_1 selection in density estimation with SIMEX}
\description{
  To compute the lambda_1 using the rule-of-thumb methods.}
\usage{
 lambda.select(W,msigma, method="SJ", na.rm = FALSE, ...)
}
\arguments{
\item{W}{The observed data.  It is a vector of length at least 10.}
\item{msigma}{The standard deviation \eqn{\sigma} of measurement error.  It is a single positive numeric value.}
\item{method}{Specifies the bandwidth of 'W'. It can be a single numeric value which has been pre-determined; or computed with the specific density bandwidth selector: 'nrd0', 'nrd', 'ucv', 'bcv', 'SJ'.}
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
 \code{\link{simex.density}}, \code{\link{simex.H.density}}.
}
\examples{

N <- 1000
set.seed(123); X <- c(rnorm(N/2, mean=-2), rnorm(N/2,mean=2)); U <- rnorm(N,sd=1)
msigma <- 0.5
#msigma <- runif(N,min=0.3,max=0.5)
W <- X + msigma*U
lambda1 <- lambda.select(W, msigma=msigma)

}
\keyword{smoothing}

