#Block Gibbs uses standard rats data 🐀
library(tidyverse)

rats <- read_csv(file="Data/rats.csv")
head(rats)
colnames(rats)

#make long
rats <- rats  %>%  gather(key="day", value="weight", day8, day15, day22, day29, day36) %>%
  mutate(id=as.factor(id),
         day=as.numeric(substring(day,4,5)))

colnames(rats)

ggplot() + geom_point(data=rats, aes(x=day, y=weight)) +
  geom_line(data=rats, aes(x=day, y=weight, group=id))

# setup random effects
y <- rats$weight
group <- rats$id
x <- rats$day


class(group)
head(group)
table(group)

Z <- model.matrix(~group-1)
head(Z)

subset(Z, group==1)



# setup fixed effects


class(x)
head(x)


# 🐀 setup design matrix
W <- cbind(1,x,Z)
head(W)

# get dimensions
q <- ncol(Z)
p <- ncol(W) - q
n <- nrow(W)

q
p
n


# setup storage
S <- 2000
beta_keep <- matrix(NA, S, p)
gamma_keep <- matrix(NA, S, q)
sig2inv_keep <- rep(NA,S)
kappa2inv_keep <- rep(NA,S)

## pre-calculations

WW <- t(W)%*%W
Wy <- t(W)%*%y

## starting values

beta <- rnorm(p) # fixed effect and intercept
gamma <- rnorm(q) # random effects only
theta <- c(beta, gamma)
sig2inv <- .1
kappa2inv <- .1



# specific covariance matrix
a1 <- a2 <- b1 <- b2 <-0.5

tau2 <- 100
SigmaDiag <- c(0, rep(1/tau2, p-1), rep(kappa2inv,q))
SigmaDiag

## MCMC

for(s in 1:S){
  ### update beta from full conditional

  # variance
  v <- WW*sig2inv
  SigmaDiag[(p+1):(p+q)] <- kappa2inv
  diag(v) <- diag(v) + SigmaDiag
  v <- chol2inv(chol(v))

  #mean
  m <- v %*% (sig2inv*Wy)

  # simulate from FC
  theta <- drop(m +  t(chol(v)) %*% rnorm(p+q))

  ### update sig2inv from full conditional
  sig2inv <- rgamma(1, a1 + n/2,  a2 + 0.5*sum((y-W%*%theta)^2))

  ### update kappa2inv from full conditional
  ### The following line was corrected after the video
  kappa2inv <- rgamma(1, b1 + q/2,  b2 + 0.5*sum(theta[(p+1):(p+q)]^2))


  ### store output
  sig2inv_keep[s] <- sig2inv
  kappa2inv_keep[s] <- kappa2inv
  beta_keep[s,] <- theta[1:p]
  gamma_keep[s,] <- theta[(p+1):(p+q)]

}

matplot(gamma_keep, type="l")

# view posterior
head(beta_keep)
matplot(beta_keep, type="l", log="y")
plot(beta_keep[,2], type="l")

sigma2 <- 1/sig2inv_keep
kappa2 <- 1/kappa2inv

# ICC

ICC <- kappa2/(kappa2+sigma2)
plot(ICC)
mean(ICC[-c(1:1000)])
quantile(ICC[-c(1:1000)], c(0.025,0.975))



mean(ICC)
mean(kappa2)/(mean(kappa2)+mean(sigma2))
