library(tidyverse)
set.seed(5106)

#simulate data
y <- rnorm(1000,5,2)
n <- length(y)
sumy <- sum(y)

# set priors
mu0 <- 0
tau2inv <- 1/10^2
alpha <- beta <- 1/2

# starting values
mu <- 0
sig2inv <- 1

R <- 1000

# place to store results
keepers <- data.frame(mu = rep(NA,R),
                   sig2inv = rep(NA,R))




# MCMC
for(s in 1:R){

  mu <- rnorm(1,
              (sig2inv*sumy + tau2inv*mu0)/(sig2inv*n + tau2inv),
              1/sqrt(sig2inv*n + tau2inv)
              )

  sig2inv <- rgamma(1,
                    alpha + n/2,
                    beta + sum((y-mu)^2)/2)

  keepers[s,] <- c(mu,sig2inv)
}

head(keepers)


# transform posterior
keepers$sigma <- 1/sqrt(keepers$sig2inv)

head(keepers)

#-----------------------------#
# diagnostics
#-----------------------------#

# trace plots
par(mfrow=c(3,1))
plot(keepers$mu, type="l", ylab="mu", xlab="MCMC iteration")
plot(keepers$sig2inv, type="l", ylab="sig2inv", xlab="MCMC iteration")
plot(keepers$sigma, type="l", ylab="sigma", xlab="MCMC iteration")

# acf plots
par(mfrow=c(3,1))
acf(keepers$mu)
acf(keepers$sig2inv)
acf(keepers$sigma)


#-----------------------------#
# summarize the posterior
#-----------------------------#


# discard first half for burn-in
keepers <- keepers[(R/2+1):R,]

# summarize posterior
mu_summary <- c(mean(keepers$mu),
  sd(keepers$mu),
  quantile(keepers$mu, 0.025),
  quantile(keepers$mu, 0.975))

sigma_summary <- c(mean(keepers$sigma),
  sd(keepers$sigma),
  quantile(keepers$sigma, 0.025),
  quantile(keepers$sigma, 0.975))

summary <- rbind(mu_summary,
      sigma_summary)

row.names(summary) <- c("mu","sigma")
colnames(summary) <- c("mean","sd","lower","upper")
round(summary,2)


#-----------------------------#
# visualize the posterior
#-----------------------------#

# plot posterior
par(mfrow=c(3,1))
hist(keepers$mu)
hist(keepers$sigma)
plot(sigma~mu, data=keepers)



ggplot(keepers, aes(x=mu)) +
  geom_density() +
  geom_vline(xintercept=quantile(keepers$mu,c(0.025,0.975))) +
  geom_vline(xintercept=mean(keepers$mu)) +
  ggtitle("Poserior distribution of mu") +
  theme_minimal()

ggplot(keepers, aes(x=sigma)) +
  geom_density() +
  geom_vline(xintercept=quantile(keepers$sigma,c(0.025,0.975))) +
  geom_vline(xintercept=mean(keepers$sigma)) +
  ggtitle("Poserior distribution of sigma") +
  theme_minimal()

ggplot(keepers, aes(x=mu, y=sigma)) +
  geom_hex()

ggplot(keepers, aes(x=mu, y=sigma)) +
  geom_point() +
  geom_density_2d()






