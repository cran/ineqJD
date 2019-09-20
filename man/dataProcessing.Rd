\name{dataProcessing}
\alias{dataProcessing}
\alias{print.dataProcessed}
\title{Data Processing}
\description{Convert raw data to frequency distribution framework and returns cumulative sums.}
\usage{
dataProcessing(
  units,
  groups = rep("G1", nrow(as.matrix(units))),
  weights = rep(1, nrow(as.matrix(units)))
)
}

\arguments{
  \item{units}{Numeric vector of length \code{n} or matrix of dimension \code{c(n, s)} containing \code{s} sources referred to \code{n} statistical units}
  \item{groups}{Vector of length \code{n} of group membership. If empty only one group is considered, otherwise the number of groups \code{g} is defined by the number of distinct values or levels in this vector.}
  \item{weights}{Vector of lweights of length \code{n}. If empty uniform weights are considered.}
}
\details{
\code{dataProcessing} convert raw data in the frequency distribution framework with \code{r} distinct values of \eqn{Y}. In this way repeated values are removed as well as ordering issues. Moreover cumulative frequencies and cumulative sources values are evaluated in order to prepare data for inequality decompositions.
}

\value{
	\item{yh}{Vector of length \code{r} of distinct values of \eqn{Y}.}
	\item{Phl}{Matrix of absolute cumulative frequencies of dimension \code{c(r, g)}.}
	\item{Qhlk}{Array of cumulative sum of sources of dimension \code{c(r, g, s)}.}
	}
\references{
Zenga M.M. (2007). Inequality curve and inequality index based on the ratios between lower and upper means. Statistica and Applicazioni, V(1), 3-27
}
\author{Alberto Arcagni, Igor Valli}


%% ~Make other sections like Warning with \section{Warning }{...} ~


\examples{
G <- c(1, 2, 3, 1, 2, 3, 1, 1, 2, 3, 3, 3) # vector denoting group membership
X1 <- c(0, 0, 0, 500, 700, 300, 750, 1000, 500, 500, 500, 1000) # vector of the first source
X2 <- c(0, 0, 0, 500, 300, 700, 750, 500, 700, 700, 1000,600) # vector of the second source
data <- data.frame(G, X1, X2) # no sample weights are considered

x <- dataProcessing(
  units = data[, c('X1', 'X2')],
  groups = data[, 'G'],
)
x
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Frequency Distribution}
