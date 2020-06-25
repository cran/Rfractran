\name{fracDo}
\alias{fracDo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The FRACTRAN Interpreter Function
}
\description{
This implements the FRACTRAN process, to wit 
For the first fraction f in the list for which nf is an integer,
 replace n by nf. 
Repeat this rule until no fraction in the list produces
 an integer when multiplied by n, then halt. 
}
\usage{
fracDo(input, fractions = NULL, nums, denoms, tries = 100,
 stopFun = NULL, liveUpdate = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
The starting integer for the given FRACTRAN code (fraction sequence) 
}
  \item{fractions}{
The sequence of fractions. Must be in \code{bigq} form. If NULL, then 
 \code{nums} and \code{denoms} must be supplied.
}
  \item{nums}{
Vector of integers (or \code{bigz} values ) representing the numerators
  of the FRACTRAN fractions. If \code{fractions} is supplied, this is ignored.
}
  \item{denoms}{
Vector of integers (or \code{bigz} values ) representing the denominators
  of the FRACTRAN fractions. If \code{fractions} is supplied, this is ignored.
}
  \item{tries}{
A 'safety' limiter on the number of times to process the next integer
 generated. This avoids possible infinite runs or other time-wasters. 
}
  \item{stopFun}{
Optional user-supplied function that can be used to terminate a FRACTRAN
 run early, or to take min-run actions such as sending information to 
  the console.  See Details for more information. 
}
  \item{liveUpdate}{
If set to TRUE, a few dots and words are sent to the console to indicate 
 that the algorithm is still running. 
}
  \item{\dots}{
Possible additional arguments for \code{stopFun} or future use. 
}
}
\details{
Some FRACTRAN programs do not terminate, most famously the prime generators. 
 If a specific value is being looked for, an appropriate \code{stopFun}
can be supplied to check for that value. \code{stopFun} must return a
 single logical value (R, as always will convert numerics to logical if
 necessary) indicating success as TRUE. The first argument
 to \code{stopFun} must accept a single \code{bigq} value.
 If there are more arguments, they must be entered after all named
 arguments.  Note that this function does not have to send TRUE; it could be 
 used solely to execute other commands mid-run. 
 See \code{\link[Rfractran]{stopPrime}} for one such example. 
 
}
\value{
A vector of all the \code{bigq} values generated. 
 These all have denominator == 1, as they are the integers found in 
 each iteration of the FRACTRAN algorithm. 
}
\references{
 \url{https://esolangs.org/wiki/Fractran}
 \url{https://oeis.org/wiki/List_of_FRACTRAN_programs_to_compute_core_sequences}
}
\author{
Carl Witthoft <carl@witthoft.com
}


\seealso{
 \code{\link[Rfractran]{fracMake}}
}
\examples{
# addition:  { 3/2   } .  enter 2^a*3^b and get 3^(a+b) answer when terminated
addit <- fracDo(2^5*3^8,nums = 3, denoms = 2)
# Last value returned is what we want. But can't take log(bigq) 
log(numerator(addit[length(addit)]),3)
}

