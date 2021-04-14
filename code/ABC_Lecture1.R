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
N <- 100
thetas <- rt(N,df=2)
unnormalizedWeights <- dnorm(thetas) / dt(thetas,df=2)
normalizedWeights <- unnormalizedWeights / sum(unnormalizedWeights)
plot(normalizedWeights)
sum(thetas * normalizedWeights)

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
    chain[s,] <- thetaStar # ACCEPT !!
  } else {
    chain[s,] <- chain[s-1,] # REJECT !!
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
chain[1,] <- 1 # good idea for a starting value

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

D <- 4
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
Y <- rnorm(N,mean = 10,sd=sqrt(2))

mu0 <- 0
tau20 <- 10

alpha <- 1
beta  <- 1

maxIts <- 1000
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
  resids <- Y-mus[i]
  sigma2s[i] <- 1/rgamma(n=1,
                         shape = alpha+N/2,
                         rate = beta + t(resids)%*%resids/2)
  
}

plot(mus,type="l")
plot(sigma2s,type="l")

#
###
#
N <- 100 # more data
Y <- rnorm(N,mean = 10,sd=sqrt(2))

mu0 <- 0
tau20 <- 10

alpha <- 1
beta  <- 1

maxIts <- 1000
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
  resids <- Y-mus[i]
  sigma2s[i] <- 1/rgamma(n=1,
                         shape = alpha+N/2,
                         rate = beta + t(resids)%*%resids/2)
  
}

plot(mus,type="l")
plot(sigma2s,type="l")  # sometimes data overpowers the prior


################################################################################
#
####
####### metropolis with adaptive covariance
####
#

recursion <- function(Ct,XbarMinus,Xt,epsilon,sd,t,warmup=100) {
  if(t>warmup) {
    XbarT  <- Xt/t + (t-1)/t*XbarMinus
    CtPlus <- Ct*(t-1)/t + sd/t*( t*XbarMinus%*%t(XbarMinus) -
                                    (t+1)*XbarT%*%t(XbarT) +
                                    Xt%*%t(Xt) + epsilon*diag(length(Xt)))
  } else {
    XbarT  <- Xt/t + (t-1)/t*XbarMinus
    CtPlus <- Ct
  }
  return(list(CtPlus,XbarT))
}

randomWalk <- function(N, x0, maxIt=10000,
                       adaptCov=FALSE) {
  if(N!=length(x0)) stop("Dimension mismatch.")
  
  chain <- matrix(0,maxIt,N)
  
  sigma <- 2.4/sqrt(N)
  Ct <- sigma^2*diag(N) 
  xbar <- x0
  
  accept <- rep(0,maxIt)
  chain[1,] <- x0
  for (i in 2:maxIt){
    if (adaptCov==FALSE) {
      xStar <- rnorm(N,sd=sigma) + chain[i-1,]
      if(log(runif(1)) < sum(target(xStar)) -
         sum(target(chain[i-1,]))){
        accept[i] <- 1
        chain[i,] <- xStar
      } else {
        chain[i,] <- chain[i-1,]
      }
    } else { # with covariance
      xStar <- as.vector(t(chol(Ct))%*%rnorm(N) + chain[i-1,])
      if(log(runif(1)) < sum(target(xStar)) -
         sum(target(as.vector(chain[i-1,])))){
        accept[i] <- 1
        chain[i,] <- xStar
      } else {
        chain[i,] <- chain[i-1,]
      }
      updt <- recursion(Ct=Ct,
                        XbarMinus=xbar,
                        Xt=chain[i,],
                        epsilon = 0.000001,
                        sd=sigma^2,
                        t=i,
                        warmup=maxIt/10)
      Ct <- updt[[1]]
      xbar <- updt[[2]]
    }
    
    if(i %% 1000 == 0) cat(i,"\n")
  }
  ratio <- sum(accept)/(maxIt-1)
  cat("Acceptance ratio: ", ratio,"\n")
  if (adaptCov) {
    return(list(chain,ratio,sigma,Ct))
  } else{
    return(list(chain,ratio,sigma,diag(N)))
  }
}

target <- function(theta) { # multivariate standard normal
  output <- mvtnorm::dmvnorm(theta, log=TRUE)
  return(output)
}

#
###
#
library(coda)


# 10 dimensions
N <- 10
results <- randomWalk(N=N,
                      x0=rep(0,N),
                      maxIt = 1000) # acceptance should be 0.234
plot(results[[1]][,1],type="l")
effectiveSize(as.mcmc(results[[1]]))

# more its
N <- 10
results <- randomWalk(N=N,
                      x0=rep(0,N),
                      maxIt = 10000) # acceptance should be 0.234
plot(results[[1]][,1],type="l")
effectiveSize(as.mcmc(results[[1]]))

# more dimensions
N <- 100
results <- randomWalk(N=N,
                      x0=rep(0,N),
                      maxIt = 10000) # acceptance should be CLOSER to 0.234
plot(results[[1]][,1],type="l")
effectiveSize(as.mcmc(results[[1]]))


#
### nasty target
#

target <- function(theta) { # multivariate standard normal
  output <- mvtnorm::dmvnorm(theta, sigma = diag(1:length(theta)), log=TRUE)
  return(output)
}

# 10 dimensions
N <- 10
results <- randomWalk(N=N,
                      x0=rep(0,N),
                      maxIt = 1000) 
plot(results[[1]][,1],type="l")
effectiveSize(as.mcmc(results[[1]]))
plot(results[[1]][,N],type="l")

# 10 dimensions
N <- 10
results <- randomWalk(N=N,
                      x0=rep(0,N),
                      maxIt = 10000) 
plot(results[[1]][,1],type="l")
effectiveSize(as.mcmc(results[[1]]))
plot(results[[1]][,N],type="l")

# 10 dimensions adapt cov
N <- 10
results <- randomWalk(N=N,
                      x0=rep(0,N),
                      adaptCov = TRUE,
                      maxIt = 10000) 
plot(results[[1]][,1],type="l")
effectiveSize(as.mcmc(results[[1]]))
plot(results[[1]][,N],type="l")

# longer chains gives more time for adaptation
N <- 10
results <- randomWalk(N=N,
                      x0=rep(0,N),
                      adaptCov = TRUE,
                      maxIt = 100000) 
plot(results[[1]][,1],type="l")
effectiveSize(as.mcmc(results[[1]]))
plot(results[[1]][,N],type="l")

# more dimensions
N <- 100
results <- randomWalk(N=N,
                      x0=rep(0,N),
                      adaptCov = FALSE,
                      maxIt = 10000) 
plot(results[[1]][,1],type="l")
effectiveSize(as.mcmc(results[[1]]))
plot(results[[1]][,N],type="l")

# now adapting
N <- 100
results <- randomWalk(N=N,
                      x0=rep(0,N),
                      adaptCov = TRUE,
                      maxIt = 10000) 
plot(results[[1]][,1],type="l")
effectiveSize(as.mcmc(results[[1]]))
plot(results[[1]][,N],type="l")

# make it beautiful
N <- 100
results <- randomWalk(N=N,
                      x0=rep(0,N),
                      adaptCov = TRUE,
                      maxIt = 100000) 
plot(results[[1]][,1],type="l")
effectiveSize(as.mcmc(results[[1]]))
plot(results[[1]][,N],type="l")


################################################################################
#
####
####### hmc for "standard" multivariate normal
####
#

library(coda)

target <- function(theta) {
  output <- - sum(theta^2)/2
  return(output)
}

grad <- function(theta) {
  output <- - theta #numDeriv::grad(target,theta)
  return(output)
}

metropolis <- function(maxIts,D) {
  chain <- matrix(0,maxIts,D) 
  chain[1,] <- 0 # good idea for a starting value
  
  for(s in 2:maxIts) {
    thetaStar <- rnorm(n=D,mean=chain[s-1,],sd=2.4/sqrt(D)) # proposal
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


hmc <- function(D, maxIts, stepSize=0.01) {

  chain <- matrix(0,maxIts,D)
  acceptances <- 0
  L <- 10 # number of leapfrog steps
  chain[1,] <- rnorm(D)
  currentU  <- - target(chain[1,])

  for (i in 2:maxIts) {
    proposalState    <- chain[i-1,]
    momentum         <- rnorm(D)
    currentK   <- sum(momentum^2)/2

    # leapfrog steps
    momentum <- momentum + 0.5 * stepSize * grad(proposalState)
    for (l in 1:L) {
      proposalState <- proposalState + stepSize * momentum
      if (l!=L) momentum <- momentum + stepSize * grad(proposalState)
    }
    momentum <- momentum + 0.5 * stepSize * grad(proposalState)
    
    # quantities for accept/reject
    proposedU = - target(proposalState)
    proposedK = sum(momentum^2)/2
    u <- runif(1)
    
    if (log(u) < currentU + currentK - proposedU + proposedK) {
      chain[i,]   <- proposalState
      currentU    <- proposedU
      acceptances <- acceptances + 1
    } else {
      chain[i,] <- chain[i-1,]
    }
    
    if (i %% 100 == 0) cat("Iteration ", i,"\n") 
  }
  
  cat("Acceptance rate: ", acceptances/(maxIts-1))
  return(chain)
}

#
### D=2
#
hmc_results <- hmc(D=2,maxIts = 10000)
rwm_results <- metropolis(D=2,maxIts = 10000)
plot(rwm_results[,1],type="l")
lines(hmc_results[,1],col="red")
effectiveSize(as.mcmc(rwm_results))
effectiveSize(as.mcmc(hmc_results)) # RWM crushes HMC!!!

#
### change hmc stepsize to 0.1
#
hmc_results <- hmc(D=2,maxIts = 10000, stepSize = 0.1)
rwm_results <- metropolis(D=2,maxIts = 10000)
plot(rwm_results[,1],type="l")
lines(hmc_results[,1],col="red")
effectiveSize(as.mcmc(rwm_results))
effectiveSize(as.mcmc(hmc_results)) # HMC crushes RWM!!!

#
### change hmc stepsize to 1
#
hmc_results <- hmc(D=2,maxIts = 10000, stepSize = 1)
rwm_results <- metropolis(D=2,maxIts = 10000)
plot(rwm_results[,1],type="l")
lines(hmc_results[,1],col="red")
effectiveSize(as.mcmc(rwm_results))
effectiveSize(as.mcmc(hmc_results)) # Wait what?!
plot(hmc_results[,1],type="l")
lines(rwm_results[,1],col="red")

#
### D=100
#
hmc_results <- hmc(D=100,maxIts = 10000, stepSize = 1)
rwm_results <- metropolis(D=100,maxIts = 10000)
effectiveSize(as.mcmc(rwm_results))
summary(effectiveSize(as.mcmc(rwm_results)))
summary(effectiveSize(as.mcmc(hmc_results))) # Wait what?!
plot(hmc_results[,1],type="l")
lines(rwm_results[,1],col="red")

#
### D=1000
#
hmc_results <- hmc(D=1000,maxIts = 10000, stepSize = 1)
rwm_results <- metropolis(D=1000,maxIts = 10000)
summary(effectiveSize(as.mcmc(rwm_results)))
summary(effectiveSize(as.mcmc(hmc_results)))
plot(hmc_results[,1],type="l")
lines(rwm_results[,1],col="red")

################################################################################
#
####
####### HMC with adaptive stepsize (0.9 acceptance rate target)
####
#

library(coda)


target <- function(theta) {
  output <- - sum(theta^2)/2
  return(output)
}

grad <- function(theta) {
  output <- - theta #numDeriv::grad(target,theta)
  return(output)
}


delta <- function(n) {
  return( min(0.01,n^(-0.5)) )
}

adapt_hmc <- function(D, maxIts, targetAccept=0.8, stepSize=1, L=20) {
  
  chain <- matrix(0,maxIts,D)
  stepSize <- 1
  chain[1,] <- rnorm(D)
  currentU  <- - target(chain[1,])
  
  totalAccept <- rep(0,maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  Proposed = 0
  
  for (i in 2:maxIts) {
    proposalState    <- chain[i-1,]
    momentum         <- rnorm(D)
    currentK   <- sum(momentum^2)/2
    
    # leapfrog steps
    momentum <- momentum + 0.5 * stepSize * grad(proposalState)
    for (l in 1:L) {
      proposalState <- proposalState + stepSize * momentum
      if (l!=L) momentum <- momentum + stepSize * grad(proposalState)
    }
    momentum <- momentum + 0.5 * stepSize * grad(proposalState)
    
    # quantities for accept/reject
    proposedU = - target(proposalState)
    proposedK = sum(momentum^2)/2
    u <- runif(1)
    
    if (log(u) < currentU + currentK - proposedU + proposedK) {
      chain[i,]   <- proposalState
      currentU    <- proposedU
      totalAccept[i] <- 1
      Acceptances = Acceptances + 1
    } else {
      chain[i,] <- chain[i-1,]
    }
    
    SampCount <- SampCount + 1

    # tune
    if (SampCount == SampBound) { 
      AcceptRatio <- Acceptances / SampBound
      if ( AcceptRatio > targetAccept ) {
        stepSize <- stepSize * (1 + delta(i-1))
      } else {
        stepSize <- stepSize * (1 - delta(i-1))
      }
  
      SampCount <- 0
      Acceptances <- 0
    }
    
    
    if (i %% 100 == 0) cat("Iteration ", i,"\n","stepSize: ", stepSize, "\n") 
  }
  
  cat("Acceptance rate: ", sum(totalAccept)/(maxIts-1))
  return(chain)
}

#
### D=1000, different target accepts
#
hmc_results  <- adapt_hmc(D=1000,maxIts = 100000,  targetAccept = 0.234)
hmc_results <- hmc_results[10001:100000,]
summary(effectiveSize(as.mcmc(hmc_results[,1:10])))
plot(hmc_results[,1],type="l")
qqnorm(hmc_results[,1])
qqline(hmc_results[,1])

hmc_results  <- adapt_hmc(D=1000,maxIts = 100000,  targetAccept = 0.9)
hmc_results <- hmc_results[10001:100000,]
summary(effectiveSize(as.mcmc(hmc_results[,1:10])))
plot(hmc_results[,1],type="l")
qqnorm(hmc_results[,1])
qqline(hmc_results[,1])

################################################################################
#
####
####### HMC with adaptive number of leapfrog steps (0.9 acceptance rate target)
####
#

library(coda)


target <- function(theta) {
  output <- - sum(theta^2)/2
  return(output)
}

grad <- function(theta) {
  output <- - theta #numDeriv::grad(target,theta)
  return(output)
}


delta <- function(n) {
  return( min(0.01,n^(-0.5)) )
}

adapt_hmc <- function(D, maxIts, targetAccept=0.8, stepSize=1, L=20) {
  
  chain <- matrix(0,maxIts,D)
  chain[1,] <- rnorm(D)
  currentU  <- - target(chain[1,])
  
  totalAccept <- rep(0,maxIts)
  Acceptances = 0 # total acceptances within adaptation run (<= SampBound)
  SampBound = 50   # current total samples before adapting radius
  SampCount = 0   # number of samples collected (adapt when = SampBound)
  Proposed = 0
  
  for (i in 2:maxIts) {
    proposalState    <- chain[i-1,]
    momentum         <- rnorm(D)
    currentK   <- sum(momentum^2)/2
    
    # leapfrog steps
    momentum <- momentum + 0.5 * stepSize * grad(proposalState)
    for (l in 1:round(L)) {
      proposalState <- proposalState + stepSize * momentum
      if (l!=round(L)) momentum <- momentum + stepSize * grad(proposalState)
    }
    momentum <- momentum + 0.5 * stepSize * grad(proposalState)
    
    # quantities for accept/reject
    proposedU = - target(proposalState)
    proposedK = sum(momentum^2)/2
    u <- runif(1)
    
    if (log(u) < currentU + currentK - proposedU + proposedK) {
      chain[i,]   <- proposalState
      currentU    <- proposedU
      totalAccept[i] <- 1
      Acceptances = Acceptances + 1
    } else {
      chain[i,] <- chain[i-1,]
    }
    
    SampCount <- SampCount + 1
    
    # tune
    if (SampCount == SampBound) { 
      AcceptRatio <- Acceptances / SampBound
      if ( AcceptRatio > targetAccept ) {
        L <- L * (1 + delta(i-1))
      } else {
        L <- L * (1 - delta(i-1))
      }
      
      SampCount <- 0
      Acceptances <- 0
    }
    
    
    if (i %% 100 == 0) cat("Iteration ", i,"\n","L: ", L, "\n") 
  }
  
  cat("Acceptance rate: ", sum(totalAccept)/(maxIts-1))
  return(chain)
}

#
### D=1000
#
hmc_results  <- adapt_hmc(D=1000,maxIts = 100000,  targetAccept = 0.9,
                          stepSize = 1.98)
hmc_results <- hmc_results[10001:100000,]
summary(effectiveSize(as.mcmc(hmc_results[,1:10])))
plot(hmc_results[,1],type="l")
qqnorm(hmc_results[,1])
qqline(hmc_results[,1])


# tune both L and stepsize?


