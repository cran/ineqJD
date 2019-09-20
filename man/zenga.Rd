\name{zenga}
\alias{zenga}
\alias{zenga.dataProcessed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Point and synthetic Zenga 2007 indexes}
\description{Computes point and synthetic Zenga 2007 indexes on a variable \eqn{Y}.}
\usage{
zenga(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x }{List containing: \code{'yh'}, the vector of unique values of the variable \eqn{Y} whose Bonferroni index is computed; \code{'Phl'}, the matrix of absoute cumulative frequencies;  \code{'Qhlk'}, the matrix of cumlative sums of \code{y} or its sources. \code{x} is usually the result of \code{dataProcessing} function. More details are given in the "Details" section and \code{\link{dataProcessing}} help page.}
}

\details{
	\code{zenga} compute point and synthetic Zenga 2007 indexes on a variable \code{y}, e.g. income, on a statistical population that could be partitioned in \code{g} subpopultions and could be considered as sum of \code{c} sources, e.g. income sources.
}
\value{
\item{index }{String denoting computed index.}
\item{decomposition }{Array containing the decompositions.}
\item{x }{Object usually of class \code{dataProcessed} passed as input.}
}
\references{

Zenga M. M.(2007). Inequality Curve and Inequality Index based on the Ratios between lower and upper means . Statistica and Applicazioni, V (1), 3-27. 

Zenga M. (2015) Joint decomposition by subpopulations and sources of the point and synthetic Zenga(2007) Index I(Y). Statistica and Applicazioni, XIII (2), pp.163-195. 
}
\author{Alberto Arcagni, Igor Valli}
%% ~Make other sections like Warning with \section{Warning }{....} ~

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
decomposition
}
\concept{Zenga concentration inequality point synthetic index decomposition}

