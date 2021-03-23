#
####
####### Sampling exponential r.v.s using quantile function
####
#
n <- 10000
lambda <- 5 # rate 

inverseCDF <- function (u,lambda) {
  return(-log(1-u)/lambda)
}

us <- runif(10000)
thetas <- inverseCDF(us,lambda)

hist(thetas)
abline(v=mean(thetas),col="red",lwd=3)

qqplot(x=thetas, y=rexp(10000,rate=5))
abline(a=0,b=1,col="red",lwd=3)


#
####
####### rejection sampler gets pi
####
#

Us <- matrix(runif(1000),500,2)
Us <- Us*2 - 1
plot(Us)

lssThn <- apply(Us, MARGIN = 1, FUN=function(x) sum(x^2)<1)
Pi <- mean(lssThn) * 4

Pis <- function (S) {
  Us <- matrix(runif(2*S),S,2)
  Us <- Us*2 - 1

  lssThn <- apply(Us, MARGIN = 1, FUN=function(x) sum(x^2)<1)
  Pi <- mean(lssThn) * 4
  
  return(Pi)
}


estimates <- unlist(lapply(5:2000, Pis))
plot(estimates,type = "l")
abline(h=pi,col="red",lwd=3)  # the more samples the better

plot(abs(estimates-pi),type = "l")
lines(sqrt(1/(5:2000)),col="red",lwd=3) # diminishing returns


#
####
####### rejection samplers for sine(0,pi)
####
#
n <- 1000
M <- 10 
proposals <- runif(n,max=pi)
envelope <- M/pi

Us <- runif(n)
acceptances <- Us < (sin(proposals)/2)/envelope
samples <- proposals[acceptances] 

plot(density(x=samples))
x <- seq(from=0,to=pi,length.out = 1000)
lines(x=x,y=sin(x)/2,col="red",lwd=3)

# on average M iterations required to generate a single sample
n/M
sum(acceptances)


rejection_sampler <- function(M,n) {
  proposals <- runif(n,max=pi)
  envelope <- M/pi
  
  Us <- runif(n)
  acceptances <- Us < (sin(proposals)/2)/envelope
  samples <- proposals[acceptances] 
  
  plot(density(x=samples),ylim=c(0,0.5))
  x <- seq(from=0,to=pi,length.out = 1000)
  lines(x=x,y=sin(x)/2,col="red",lwd=3)
  
  cat("n/M:", n/M, ", Number of samples: ",sum(acceptances),"\n")
}


rejection_sampler(n=1000,M=100)
rejection_sampler(n=1000,M=10)
rejection_sampler(n=1000,M=5)
rejection_sampler(n=1000,M=pi) 
rejection_sampler(n=1000,M=pi/2) # choose a good envelope


rejection_sampler(n=10000,M=pi/2)
rejection_sampler(n=100000,M=pi/2)
rejection_sampler(n=1000000,M=pi/2) # and collect many samples



#
####
####### importance sampler: Gaussian and t-distribution
####
#

