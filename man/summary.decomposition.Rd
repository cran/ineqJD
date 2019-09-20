\name{summary.decomposition}
\alias{summary.decomposition}
\alias{print.summary.decomposition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing inequality decomposition}

\description{
  \code{summary} method for class \code{"decomposition"}.
}
\usage{
\method{summary}{decomposition}(object, ...)
\method{print}{summary.decomposition}(x, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{"decomposition"}, usually, as result of a call to \code{\link{gini}}, \code{\link{bonferroni}} and \code{\link{zenga}}.
%%     ~~Describe \code{x} here~~
}
  \item{x}{rtrtrt}
  \item{\dots}{further arguments passed to or from other methods.}
%%     ~~Describe \code{\dots} here~~
}
\details{
\code{summary.decomposition} method use 
}
\value{

  \item{index }{String denoting computed index.}
  \item{joint }{Array of joint decompositions by sources and subpopulations.}
  \item{pairs }{Matrix of decompositions by subpopulations.}
  \item{within }{Vector of \code{within} part to the overall inequality. It denotes the part of the overall inequality derived from the inequality inside each subpopulation.}
  \item{between }{Vector of \code{between} part to the overall inequality. It denotes the part of the overall inequality derived from the comparison between subpopulations.}
  \item{groups }{Vector of subpopulations contribution to the overall inequality.}
  \item{groups_sources }{Matrix of subpopulations contributions for each source to the overall inequality.}
  \item{sources }{Vector of sources contribution to the overall inequality.}
  \item{synthetic }{Scalar denoting the value of the synthetic index.}
%% ...
}
\references{

Zenga M. M.(2007). Inequality Curve and Inequality Index based on the Ratios between llower and upper Means . Statistica and Applicazioni, V (1), 3-27. 

Zenga M. (2015) Joint decomposition by subpopulations and sources of the point and synthetic Zenga(2007) Index I(Y). Statistica and Applicazioni, XIII (2), pp.163-195.

Zenga M., Valli I. (2017). Joint decomposition by Subpopulations and Sources of the Point and Synthetic Bonferroni Inequality Measures. Statistics and Applications, XV (2), pp. 83-120.

Zenga M., Valli I. (2018). Joint decomposition by Subpopulations and Sources of the Point and Synthetic Gini Indexes. Statistics and Applications, XVI (1).

}
\author{
Alberto Arcagni, Igor Valli.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{gini}}, \code{\link{bonferroni}}, \code{\link{zenga}}, \code{\link{dataProcessing}}.
}
\examples{
G <- c(1, 2, 3, 1, 2, 3, 1, 1, 2, 3, 3, 3) # vector denoting group membership
X1 <- c(0, 0, 0, 500, 700, 300, 750, 1000, 500, 500, 500, 1000) # vector of the first source
X2 <- c(0, 0, 0, 500, 300, 700, 750, 500, 700, 700, 1000,600) # vector of the second source
data <- data.frame(G, X1, X2) # no sample weights are considered
x <- dataProcessing( # data preparation
  units = data[, c('X1', 'X2')],
  groups = data[, 'G'],
)
decomposition <- zenga(x)

summary(decomposition)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

