\name{L1.test}
\alias{L1.test}
\title{
  L1 permutation test
}
\description{
  Two-sample permutation test with two-sided alternative with the following statistics:
  \deqn{\sum[1\le i \le n.x]\sum[1\le j \le n.y] ln(1+|x[i]-y[j]|)}
}
\usage{
  L1.test(x, y, normalization = FALSE, n.permutations = 10000, exact = NULL, permutations = NULL, check_permutations = TRUE, ...)
}
\arguments{
  \item{x}{
    a (non-empty) numeric vector of data values.
  }
  \item{y}{
    a (non-empty) numeric vector of data values.
  }
  \item{normalization}{
    a logical indicating whether to use samples scaling. The normalizator is
    \deqn{2\sum[1\le i<j\le n] |z[i]-z[j]|/(n(n-1)), z = c(x,y), n=n.x+n.y.}
  }
  \item{n.permutations}{
    a number indicating the number of random permutations. Used only if "permutations" argument is not given.
  }
  \item{exact}{
    a logical indicating whether to use a permutation set equal to the entire permutation set. If NULL automatically set to TRUE if power of exact permutation set less or equal  "n.permutations" argument. Used only if "permutations" argument is not given.
  }
  \item{permutations}{
    an optional matrix containing permutation set. Permutations should place in the columns of the matrix and contain values coincided with the given samples.
  }
  \item{check_permutations}{
    an optional logical indicating whether to check the given permutation set. Used only if "permutations" argument is given.
  }
  \item{\dots}{
    additional parameters.
  }
}
\value{
  A list with class "htest" containing the following components:
  \item{p.value}{the p-value for the test.}
  \item{normalization}{a logical indication of whether normalization has been used.}
  \item{normalizer}{value of used normalizer.}
  \item{permutations}{permutation set that has been used.}
  \item{n.permutations}{size of used permutation set.}
  \item{exact}{a logical indication of whether exact permutation set has been used.}
  \item{n.x}{length of x sample}
  \item{n.y}{length of y sample}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string indicating what type of t-test was performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\author{
  Dmitry Salnikov (mejibkop.ru@mail.ru)
}
\examples{
  x <- rnorm(5)
  y <- rnorm(10)
  
  L1.test(x, y)
}