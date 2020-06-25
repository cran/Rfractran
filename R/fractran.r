# for the first fraction f in the list for which nf is an integer, replace n by nf
# repeat this rule until no fraction in the list produces an integer when multiplied by n, then halt.

# validation tests:
# addition:  { 3/2   } .  enter 2^a*3^b and get 3^(a+b) answer when terminated

# multiplication: complicated because the sequence "bounces" between two subsets of the fraction set. Notice that the bigger denom, 33=11*3 in this case, ALWAYS has to be listed first before the denom 3, or it will never be reached. 
# { 455/33, 11/13, 1/11, ; 3/7, 11/2, 1/3 } ( ";" there for emphasis only)

# then  2^a*3^b  terminates with  5^(a*b) 
# fracmul <-makefrac(c(455,11,1,3,11,1),c(33,13,11,7,2,3))


# prime generator; moves a lot faster/farther than Conway's 
#start with 10, feed {7/3 99/98 13/49 39/35 36/91 10/143 49/13 7/11 1/2 91/1} and result includes a series of values 10^p , p prime 
# fracp10 <- fracMake(c(7,99,13,39,36,10,49,7,1,91), c(3,98,49,35,91,143,13,11,2,1)) 

#Conway's original prime gen : feed it 2, get  2^(primes) here and ther
#fracpc <-makefrac(c(17,78,19,23,29,77,95,77,1,11,13,15,1,55), c(91,85,51,38,33,29,23,19,17,13,11,2,7,1))

#Remember - if you run N 'tries' you can restart with the Nth output integer. 




# first, a "set maker" which takes  n,d vectors or matrix and makes the right bigq out of them .
fracMake <- function(n, d = NULL){
#notice that entering a vector 'n' and no 'd' will produce the sequence 1/n[1],2/n[2], ...
nd <-xy.coords(as.bigz(n),as.bigz(d))
return(invisible(as.bigq(nd$x,nd$y)))
}

# A Q&D  logulator to pull numbers of interest out of the returned set of bigqs 

fracAns <- function(intvec, logbase, logprec = 1e-5) {
logall <- log(numerator(intvec),logbase)
frist <- which(abs(logall - round(logall)) < logprec )  # use 1e-4 as a solid place
return(invisible(logall[frist]))

}

#	2) Figure out how to get it at least 10X faster.
#	1) Notice that this returns nothing (empty) if the value of 'input' is a winner.
#		Change so first returned value is the input value 
# Here's the fraccinator engine.
# 'tries' is default terminator so things don't go blooey 
fracDo <-function(input, fractions=NULL, nums, denoms, tries = 100, stopFun=NULL, liveUpdate=TRUE , ...) {
input <- as.bigq(input) # for safety, do everything in bigq rather than bigz
if (length(fractions)) {
set <- fractions
if (!is.bigq(set)) stop('fraction set must be bigq')
} else {
	set <-as.bigq(floor(nums), floor(denoms))
}
setlen <- length(set)
intvec = input # will store 'wins' here # but bigq(NULL) is weird.  
itry = 1
while (itry <= tries) {
	if (liveUpdate) {
		if(!itry%%1000){
			cat(' .') 
# moved inside to reduce if-loading
			if (!itry%%20000) {
				cat('\n')
#				flush.console()
			}
			flush.console()
		}  #, appendLF = FALSE)
	}
	idx = 1 # index for grabbing a fraction 
	gotit = 1 # set to zero if new integer found
	while (gotit) {
		tmpterm <- mul.bigq(input, set[idx])
		if (denominator(tmpterm ) == 1){
			intvec <- c(intvec,tmpterm)
# if termination desired, do it here
			if(!is.null(stopFun)) {
				foundit <- stopFun(tmpterm, ...) #must return a single logical
				flush.console()
				if (foundit){
					return(invisible(intvec))
				}
			}
			gotit <- 0
			input <- tmpterm # start again with new integer
		} else {
			idx <- idx + 1
		}
		if (idx > setlen){
			gotit = 0
			itry <- tries + 1 # terminate because didn't find another int
		}
	} #end of gotit
	itry = itry + 1 # work on new "input" integer
} # end of tries
if (liveUpdate) cat('\n') # just to get a newline
return(invisible(intvec ))
}

# Let's create a "stopFun" whose only job is to cat() a found prime without stopping the main function.
# TODO: 
#2: see if simple unique(factorize) test works better than my goofball code 
#1:  modify to work with either the 10^prime or the 2^prime code.
# Logb x = Loga x/Loga b  so log(x,mant) = log(x)/log(mant)
# this will require 'mantissa' to be in the '...' of fracDo .See if neater way

# cheap time test.  catPrime:     
# user  system elapsed 
#  64.613   0.712  65.737 
# stopPrime
#  user  system elapsed 
#  63.588   0.398  63.879 
 
stopPrime <- function(x,mantissa,returnVal = 0) {
mantissa = round(mantissa) # safety
manfac <- factorize(mantissa)
facs <- factorize(numerator(x)) # x is a bigq
# cyclically remove sets of 'facs' which make up one manfac
if (prod(unique(manfac) == unique(facs)) ) {
# no extraneous values in factorization
	manle <-rle(as.numeric(manfac))
	numle <- rle(as.numeric(facs))
# see if all factors are proportional quantities
	if (length(unique(manle$lengths/numle$lengths)) == 1) {
		 cat(c('found ',log(numerator(x) )/log(mantissa),'\n') )
		flush.console()
	} 
}
return(invisible(returnVal)) # typically don't want to terminate fracDo
}

# catPrime <-function(x){
# # x will always be a bigq
# aa <- as.numeric(unlist(strsplit(as.character(numerator(x)),'')))
# if (aa[1] == 1 && sum(aa[2:length(aa)]) == 0) {
# 	 cat(c('found ',length(aa)-1),'\n')
# 	flush.console()
# 	} 
# return(invisible(0)) # Since I don't want to terminate fracDo
# }

