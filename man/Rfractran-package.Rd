\name{Rfractran-package}
\alias{Rfractran-package}
\alias{Rfractran}
\docType{package}
\title{
\packageTitle{Rfractran}
}
\description{
\packageDescription{Rfractran}
}
\details{
The DESCRIPTION file:
\packageDESCRIPTION{Rfractran}
\packageIndices{Rfractran}
FRACTRAN is an obscure yet tantalizing programming "language"
 invented by John Conway of "Game of Life" fame. The code consists
 of a sequence of fractions.  The operation is simple.
 1 - Initalize with an integer
 2 - Multiply the integer by the first fraction. If an integer results,
 start again with the new integer. If not, try the next fraction.
 3 - If no such multiplication yields an integer, terminate the program.
 One warning: there is a FRACTRAN program that can be found on the web
  which is supposed to generate the digits of pi.  Unfortunately, it's 
  known to have a bug which causes it not to work. 
 So far nobody has found a correction for it.
}
\author{
\packageAuthor{Rfractran}
Maintainer: \packageMaintainer{Rfractran}
}
\references{
 \url{https://esolangs.org/wiki/Fractran}
 \url{https://oeis.org/wiki/List_of_FRACTRAN_programs_to_compute_core_sequences}
}
\examples{
##Not Run
# addition:  { 3/2 } . Enter 2^a*3^b and get 3^(a+b) answer when terminated

# multiplication:  
# { 455/33, 11/13, 1/11,  3/7, 11/2, 1/3 } 
# then  enter 2^a*3^b  Terminates with  5^(a*b)

# prime generator (This function never terminates.)
#start with 10, feed to [7/3 99/98 13/49 39/35 36/91 10/143 49/13 7/11 1/2 91/1]
# and whenever the result of a run is of the form 10^p , p is prime 
}