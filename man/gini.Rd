\name{gini}
\alias{gini}
\alias{gini.dataProcessed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Point and synthetic Gini indexes}
\description{Computes point and synthetic Gini indexes on a variable \eqn{Y}.}
\usage{
gini(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x }{List containing: \code{'yh'}, the vector of unique values of the variable \eqn{Y}   whose Bonferroni index is computed; \code{'Phl'}, the matrix of absoute cumulative frequencies;  \code{'Qhlk'}, the matrix of cumlative sums of \eqn{Y} or its sources. \code{x} is usually the result of \code{dataProcessing} function. More details are given in the "Details" section and \code{\link{dataProcessing}} help page.}
}
\details{
	\code{gini} compute point and synthetic Gini indexes on a variable \code{y}, e.g. income, on a statistical population that could be partitioned in \code{g} subpopultions and could be considered as sum of \code{c} sources, e.g. income sources.
}
\value{
\item{index }{String denoting computed index.}
\item{decomposition }{array containing the decompositions.}
\item{x }{object usually of class \code{dataProcessed} passed as input.}
}
\references{
Zenga M., Valli I. (2018). Joint decomposition by Subpopulations and Sources of the Point and Synthetic Gini Indexes. Statistics and Applications, XVI (1).
}
\author{Alberto Arcagni, Igor Valli.}


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
  
decomposition <- gini(x)
decomposition
}
\concept{Gini concentration inequality point synthetic index decomposition}