\name{stopPrime}
\alias{stopPrime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to Stop \code{fracDo} When Finding Primes
}
\description{
This function is designed to be fed to \code{fracDo} as the argument
 \code{stopFun}. It will send a note to the console every time a prime number
 is found. It also serves as a sample of what a "stopping function" could 
 look like - notice that despite the name, this function doesn't 
 force-stop \code{fracDo}.
}
\usage{
stopPrime(x, mantissa,returnVal = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The current integer generated with the FRACTRAN generator (\code{fracDo}.)
}
  \item{mantissa}{
Choose the mantissa based on the prime-generating algorithm selected.
Must be an added argument to \code{fracDo}. See Details. 
}
\item{returnVal}{Normally set to zero so as not to terminate
the FRACTRAN program. If desired, set to one (or TRUE) to terminate
after the first prime is found. Must be an added argument to \code{fracDo}. See Details. 
}
}
\details{
There are two well-known prime generators for FRACTRAN, one of which uses 
 powers of two and the other uses powers of 10. Select the mantissa to match 
 one of these, or potentially any other FRACTRAN program which returns 
 values as (mantissa)^(result) .  When running \code{fracDo}, you must
 enter the value as an argument after all named arguments. Similarly, if
 desired, enter a value for \code{returnVal} to \code{fracDo}.
 }
\value{
Returns the specified value of \code{returnVal} so that \code{fracDo}
 will or will not terminate. Its "useful" output is cat-ed to the console.
}

\author{
Carl Witthoft <carl@witthoft.com>
}

\examples{

##---The prime generator doesn't terminate, so look for values.
#start with 10,  and result includes a series of values 10^p , p prime 
fracp10 <- fracMake(c(7,99,13,39,36,10,49,7,1,91), c(3,98,49,35,91,143,13,11,2,1)) 
p10 <-fracDo(10,fractions = fracp10, tries = 1000,stopFun=stopPrime, mantissa = 10 )
# primes found are sent to console 
# p10 contains all integers found.
}