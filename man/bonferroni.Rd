\name{bonferroni}
\alias{bonferroni}
\alias{bonferroni.dataProcessed}

\title{Point and synthetic Bonferroni indexes}
\description{Computes the decomposition of the Bonferroni point inequality indexes of a statistical variable \eqn{Y} described in the object \code{x}.}
\usage{
bonferroni(x)
}

\arguments{
  \item{x }{An object of class \code{"dataProcessed"}. \code{x} is usually the result of \code{dataProcessing} function. More details are given in the "Details" section and \code{\link{dataProcessing}} help page.}
}

\details{
	\code{bonferroni} computes the decomposition of the Bonferroni point inequality indexes from the object \code{x} of class \code{"dataProcessed"}. \code{x} is usually the result of \code{dataProcessing} function.
}
\value{
\item{index }{String denoting computed index.}
\item{decomposition }{Array containing the decompositions. The dimensions of \code{decomposition} are \code{c(g, g, r, s)} where \code{g} is the number of groups, \code{r} the number of different values of \eqn{Y} and \code{s} the number of sources.}
\item{x }{Object of class \code{dataProcessed} passed as input.}
}
\references{
	Zenga M., Valli I. (2017). Joint decomposition by Subpopulations and Sources of the Point and Synthetic Bonferroni Inequality Measures. Statistics and Applications, XV (2), pp. 83-120. 
}
\author{Alberto Arcagni, Igor Valli.}
\seealso{
\code{\link{gini}} and \code{\link{zenga}} for other inequality indexes and \code{\link{dataProcessing}} for the class \code{"dataProcessed"}.
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
  
decomposition <- bonferroni(x)
decomposition
}
\concept{Bonferroni concentration inequality point synthetic index decomposition}