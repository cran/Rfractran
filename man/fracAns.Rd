\name{fracAns}
\alias{fracAns}

\title{
A Q&D  "Logulator" to Find Numbers of Interest From \code{fracDo} Output
}
\description{
In nearly all cases FRACTRAN's "result" is the exponent of some number.
  This function takes the logarithm of the specified base and 
  identifies values which are integers (or nearly so, to the specified precision).
}
\usage{
fracAns(intvec, logbase, logprec = 1e-05)
}

\arguments{
  \item{intvec}{
A vector of \code{bigq} values, as returned from \code{fracDo}
}
  \item{logbase}{
The base of the desired logarithm, e.g. 2 or 10 in many cases. 
}
  \item{logprec}{
A reasonably small value used to check whether the returned double should be considered to be an integer (thus ignoring binary precision errors)
}
}

\value{
A vector of the integer values found
}

\author{
Carl Witthoft  <carl@witthoft.com>
}

\examples{
##---The prime generator doesn't terminate, so look for values.
#start with 10,  and result includes a series of values 10^p , p prime 
fracp10 <- fracMake(c(7,99,13,39,36,10,49,7,1,91), c(3,98,49,35,91,143,13,11,2,1)) 
p10 <-fracDo(10,fractions = fracp10, tries = 1000)
foundp10 <-fracAns(p10,logbase = 10)
# [1] 1 2 3 5 7 # sorry about the "1" :-) 

}