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

################################################################################
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

################################################################################
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


################################################################################
#
####
####### importance sampler: Gaussian and t-distribution
####
#

# first use gaussian to sample t-distribution
N <- 100
thetas <- rnorm(N)
unnormalizedWeights <-  dt(thetas,df=2)/dnorm(thetas)
normalizedWeights <- unnormalizedWeights / sum(unnormalizedWeights)
plot(normalizedWeights)

estimate <- sum(thetas * normalizedWeights) # compute mean of t distribution (0)
estimate

# make function and get variance of estimator
gaussiansTargetT <- function(N,maxIts) {
  estimates <- rep(0,maxIts)
  for(s in 1:maxIts) {
    thetas <- rnorm(N)
    unnormalizedWeights <-  dt(thetas,df=2)/dnorm(thetas)
    normalizedWeights <- unnormalizedWeights / sum(unnormalizedWeights)
    estimates[s] <- sum(thetas * normalizedWeights) # compute mean of t distribution (0)
  }

  return(estimates)
}
ests <- gaussiansTargetT(N=100,maxIts=1000)
plot(density(ests))
var(ests)

# now t targets gaussian
tTargetsGaussian <- function(N,maxIts) {
  estimates <- rep(0,maxIts)
  for(s in 1:maxIts) {
    thetas <- rt(N,df=2)
    unnormalizedWeights <- dnorm(thetas) / dt(thetas,df=2)
    normalizedWeights <- unnormalizedWeights / sum(unnormalizedWeights)
    estimates[s] <- sum(thetas * normalizedWeights) # compute mean of t distribution (0)
  }
  
  return(estimates)
}
ests <- tTargetsGaussian(N=100,maxIts=1000)
plot(density(ests))
var(ests)


################################################################################
#
####
####### curse of dimensionality (importance sampling)
####
#

library(mvtnorm)

# high dimensional probability density functions are small
x <- rnorm(100)
pdfs <- dnorm(x)
joint_pdf <- prod(pdfs)

x <- rnorm(600)
pdfs <- dnorm(x)
joint_pdf <- prod(pdfs)

# ratios of small numbers do not behave well
thetas <- rmvnorm(n=1000,sigma = diag(100),mean = rep(1,100))
dim(thetas)

unnormalizedWeights <-  dmvnorm(thetas)/dmvnorm(thetas,mean = rep(1,100))
normalizedWeights <- unnormalizedWeights / sum(unnormalizedWeights)
plot(normalizedWeights)
plot(density(normalizedWeights))
summary(normalizedWeights) # underflow

estimate <- colSums(thetas * normalizedWeights) # compute mean of t normals (0)
plot(density(estimate))

# higher dim = 200
thetas <- rmvnorm(n=1000,sigma = diag(200),mean = rep(1,200))
dim(thetas)

unnormalizedWeights <-  dmvnorm(thetas)/dmvnorm(thetas,mean = rep(1,200))
normalizedWeights <- unnormalizedWeights / sum(unnormalizedWeights)
plot(normalizedWeights)
plot(density(normalizedWeights))
summary(normalizedWeights) # underflow

estimate <- colSums(thetas * normalizedWeights) # compute mean of t normals (0)
plot(density(estimate))

# higher dim = 200 / overpower with 10000 samples
thetas <- rmvnorm(n=10000,sigma = diag(200),mean = rep(1,200))
dim(thetas)

unnormalizedWeights <-  dmvnorm(thetas)/dmvnorm(thetas,mean = rep(1,200))
normalizedWeights <- unnormalizedWeights / sum(unnormalizedWeights)
plot(normalizedWeights)
plot(density(normalizedWeights))
summary(normalizedWeights) # underflow

estimate <- colSums(thetas * normalizedWeights) # compute mean of t normals (0)
plot(density(estimate))

################################################################################
#
####
####### finite markov chain
####
#
library(markovchain)

Q <- matrix(c(0,0,0.5,1,0.2,0.5,0,0.8,0),3,3)  # transition matrix
mc <- new('markovchain',
                     transitionMatrix = Q, 
                     states = c('A','B','C'))

layout <- matrix(c(-2,0,0,1,2,0), ncol = 2, byrow = TRUE)
plot(mc, layout = layout)

# rows sum to 1
Q %*% c(1,1,1)

# marginal probability after 1000 steps
p <- c(1,1,1)/3
for (i in 1:1000) {
  p <- p %*% Q
}
p
sum(p) # sanity check
p == p %*% Q # is stationary distribution?

# different starting values
p <- c(1,0,0)
for (i in 1:1000) {
  p <- p %*% Q
}
p
sum(p) # sanity check
p == p %*% Q 

eig.obj <- eigen(t(Q))
eig.obj$values # hard to read
Mod(eig.obj$values)

- eig.obj$vectors * sqrt(sum(p^2)) # r scales to norm 1


# actually simulate from the markov chain 
n <- 10000
chain <- rep(0,n)
chain[1] <- 1
for (i in 2:n) {
  chain[i] <- sample(1:3,size = 1,prob = Q[chain[i-1],])
}
plot(chain) # not instructive
hist(chain)
p*n # stationary probabilities * iterations

# get mean of f(x)=x^2
1/n*sum(chain^2)

# get true mean
sum(p*c(1,2,3)^2)

# so, what's the point?

################################################################################
#
####
####### metropolis
####
#

D <- 4
maxIts <- 1000

# target function is proportional to log pdf of D-dimensional spherical Gaussian
target <- function(theta) {
  output <- mvtnorm::dmvnorm(theta, log=TRUE)
  return(output)
}

chain <- matrix(0,maxIts,D) 
chain[1,] <- - 10 # bad idea for a starting value

for(s in 2:maxIts) {
  thetaStar <- rnorm(n=D,mean=chain[s-1,])
  u         <- runif(1)
  
  logA      <- target(thetaStar) - target(chain[s-1,]) # target on log scale
  
  if(log(u) < logA) {
    chain[s,] <- thetaStar
  } else {
    chain[s,] <- chain[s-1,]
  }
}

plot(chain[,1],type="l") # not bad mixing
plot(density(chain[,1])) # meh

# make function and try with different maxIts

metropolis <- function(maxIts,D) {
  chain <- matrix(0,maxIts,D) 
  chain[1,] <- - 10 # bad idea for a starting value
  
  for(s in 2:maxIts) {
    thetaStar <- rnorm(n=D,mean=chain[s-1,]) # proposal
    u         <- runif(1)
    
    logA      <- target(thetaStar) - target(chain[s-1,]) # target on log scale
    if(log(u) < logA) {
      chain[s,] <- thetaStar
    } else {
      chain[s,] <- chain[s-1,]
    }
    
    if(s %% 100 == 0) cat(s,"\n")
  }
  
  return(chain)
}

results <- metropolis(10000, D)
plot(results[,1],type="l")
plot(density(results[,1]))


results <- metropolis(100000, D)
plot(results[,1],type="l")
plot(density(results[,1]))

D <- 100 # curse of dimensionality
results <- metropolis(10000, D)
plot(results[,1],type="l")
plot(density(results[,1]))  # how to choose proposal variance?


################################################################################
#
####
####### metropolis hastings (for truncated normal target)
####
#

D <- 4
maxIts <- 1000

# target function is proportional to log pdf of D-dimensional spherical Gaussian
target <- function(theta) {
  output <- mvtnorm::dmvnorm(theta, log=TRUE)
  return(output)
}

chain <- matrix(0,maxIts,D) 
chain[1,] <- 50 # bad idea for a starting value

for(s in 2:maxIts) {
  thetaStar <- truncnorm::rtruncnorm(n=D, a=0, b=Inf, mean = chain[s-1,], sd = 1)
  u         <- runif(1)
  
  logA      <- target(thetaStar) - target(chain[s-1,]) + # targets
               sum(log(truncnorm::dtruncnorm(x=chain[s-1,],a=0, mean=thetaStar))) -
               sum(log(truncnorm::dtruncnorm(x=thetaStar,a=0, mean=chain[s-1,])))
                
  if(log(u) < logA) {
    chain[s,] <- thetaStar
  } else {
    chain[s,] <- chain[s-1,]
  }
}

plot(chain[,1],type="l") # finds the high density area
hist(chain[,1],freq = FALSE,breaks = 20)
x <- seq(from=0,to=4,length.out = 1000)
lines(x=x,y=truncnorm::dtruncnorm(x,a=0),col="red",lwd=3)

# now with better starting values
chain <- matrix(0,maxIts,D) 
chain[1,] <- 1 # bad idea for a starting value

for(s in 2:maxIts) {
  thetaStar <- truncnorm::rtruncnorm(n=D, a=0, b=Inf, mean = chain[s-1,], sd = 1)
  u         <- runif(1)
  
  logA      <- target(thetaStar) - target(chain[s-1,]) + # targets
    sum(log(truncnorm::dtruncnorm(x=chain[s-1,],a=0, mean=thetaStar))) -
    sum(log(truncnorm::dtruncnorm(x=thetaStar,a=0, mean=chain[s-1,])))
  
  if(log(u) < logA) {
    chain[s,] <- thetaStar
  } else {
    chain[s,] <- chain[s-1,]
  }
}

plot(chain[,1],type="l") # ok!
hist(chain[,1],freq = FALSE,breaks = 20)
lines(x=x,y=truncnorm::dtruncnorm(x,a=0),col="red",lwd=3)

# make function and try with different maxIts

metropolis_hastings <- function(maxIts,D) {
  chain <- matrix(0,maxIts,D) 
  chain[1,] <- 1 
  
  for(s in 2:maxIts) {
    thetaStar <- truncnorm::rtruncnorm(n=D, a=0, b=Inf, mean = chain[s-1,], sd = 1)
    u         <- runif(1)
    
    logA      <- target(thetaStar) - target(chain[s-1,]) + # targets
      sum(log(truncnorm::dtruncnorm(x=chain[s-1,],a=0, mean=thetaStar))) -
      sum(log(truncnorm::dtruncnorm(x=thetaStar,a=0, mean=chain[s-1,])))
    
    if(log(u) < logA) {
      chain[s,] <- thetaStar
    } else {
      chain[s,] <- chain[s-1,]
    }
    
    if(s %% 100 == 0) cat(s,"\n")
  }
  
  return(chain)
}

results <- metropolis_hastings(10000, D)
plot(results[,1],type="l")
hist(results[,1],freq = FALSE,breaks = 20)
lines(x=x,y=truncnorm::dtruncnorm(x,a=0),col="red",lwd=3)

results <- metropolis_hastings(100000, D)
plot(results[,1],type="l")
hist(results[,1],freq = FALSE,breaks = 20)
lines(x=x,y=truncnorm::dtruncnorm(x,a=0),col="red",lwd=3)

D <- 100 # curse of dimensionality
results <- metropolis_hastings(10000, 20)
hist(results[,1],freq = FALSE,breaks = 20)
lines(x=x,y=truncnorm::dtruncnorm(x,a=0),col="red",lwd=3)

################################################################################
#
####
####### gibbs sampler for simplest Gaussian model
####
#

N <- 10 # data
Y <- rnorm(N,mean = 10)

mu0 <- 0
tau20 <- 100

alpha <- 1
beta  <- 1

maxIts <- 100
mus    <- rep(0,maxIts) # chains
mus[1] <- mu0
sigma2s <- rep(0,maxIts)
sigma2s[1] <- tau20

for (i in 2:maxIts) {
  # mu update
  Var <- 1/(1/tau20 + N/sigma2s[i-1])
  Mean <- (mu0/tau20+ sum(Y)/sigma2s[i-1])*Var
  mus[i] <- rnorm(n=1,mean=Mean,sd=sqrt(Var))
  
  # sigma2 update
  
  
}

