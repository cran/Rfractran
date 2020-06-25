\name{fracMake}
\alias{fracMake}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to Create a Sequence of \code{bigq} Fractions
}
\description{
Feed this function a collection of numerator and denominator values;
 get a vector of \code{bigq} fractions suitable for use as a FRACTRAN program.
 }
\usage{
fracMake(n, d = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Vector of integers to be used as numerators, or a Nx2 array of integers
 where each row contains a num, denom pair. If this is an arry, the input
  \code{d} is ignored.
}
  \item{d}{
Vector of integers to be used as denominators. Ignored if \code{n} is an array
}
}

\value{
Vector of \code{bigq} fractions. 
}

\author{
Carl Witthoft <carl@witthoft.com>
}
