library(tidyverse)
library(jagshelper)

# STILL NEED TO DEAL WITH DATES FORMATTED AS "9-May" etc

## maybe the better solution is to pivot_longer and ggplot

sonar <- read.csv("flat_data/Copper River sockeye sonar counts 1980-2025.csv")
largefish <- read.csv("flat_data/LargeFish1.csv")
temperature <- read.csv("flat_data/Temperature0.csv") %>%
  filter(`Useable.`==1)


## Investigating sonar (ALL FISH)
sonar_cumul <-
  sonar_cumulprop <-
  sonar0 <-
  sonar # initializing
for(j in 2:ncol(sonar)) {
  sonar0[,j] <- ifelse(is.na(sonar[,j]), 0, sonar[,j])
  sonar_cumul[,j] <- cumsum(sonar0[,j])
  sonar_cumulprop[,j] <- sonar_cumul[,j]/sum(sonar[,j], na.rm=TRUE)
}

cols <- rcolors(100)
# plotlines <- function(x, col=NA) {
#   if(is.na(col)) col <- rcolors(ncol(x))
#   plot()
# }

plot(NA,
     ylim=range(0, sonar[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(sonar)))
for(j in 2:ncol(sonar)) lines(sonar[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(sonar[,-1], na.rm=TRUE), lwd=3)

plot(NA,
     ylim=range(0, sonar_cumul[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(sonar_cumul)))
for(j in 2:ncol(sonar_cumul)) lines(sonar_cumul[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(sonar_cumul[,-1], na.rm=TRUE), lwd=3)

plot(NA,
     ylim=range(0, sonar_cumulprop[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(sonar_cumulprop)))
for(j in 2:ncol(sonar_cumulprop)) lines(sonar_cumulprop[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(sonar_cumulprop[,-1], na.rm=TRUE), lwd=3)




## Investigating LARGE FISH
largefish_cumul <-
  largefish_cumulprop <-
  largefish0 <-
  largefish # initializing
for(j in 2:ncol(largefish)) {
  largefish0[,j] <- ifelse(is.na(largefish[,j]), 0,
                           ifelse(largefish[,j] < 0, 0,
                                  largefish[,j]))
  largefish_cumul[,j] <- cumsum(largefish0[,j])
  largefish_cumulprop[,j] <- largefish_cumul[,j]/sum(largefish[,j], na.rm=TRUE)
}

cols <- rcolors(100)
# plotlines <- function(x, col=NA) {
#   if(is.na(col)) col <- rcolors(ncol(x))
#   plot()
# }

plot(NA,
     ylim=range(0, largefish[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(largefish)))
for(j in 2:ncol(largefish)) lines(largefish[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(largefish[,-1], na.rm=TRUE), lwd=3)

plot(NA,
     ylim=range(0, largefish_cumul[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(largefish_cumul)))
for(j in 2:ncol(largefish_cumul)) lines(largefish_cumul[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(largefish_cumul[,-1], na.rm=TRUE), lwd=3)

plot(NA,
     ylim=range(0, largefish_cumulprop[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(largefish_cumulprop)))
for(j in 2:ncol(largefish_cumulprop)) lines(largefish_cumulprop[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(largefish_cumulprop[,-1], na.rm=TRUE), lwd=3)



## Investigating LARGE FISH PROPORTION
largefish_prop <- largefish #initializing
for(j in 2:ncol(largefish)) {
  largefish_prop[,j] <- largefish[,j]/sonar[[names(largefish)[j]]]
}
plot(NA,
     ylim=range(0, largefish_prop[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(largefish_prop)))
for(j in 2:ncol(largefish_prop)) lines(largefish_prop[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(largefish_prop[,-1], na.rm=FALSE), lwd=3)

## try a hierarchical logistic regression thingy
## - maybe this would actually be robust to lag??





# do run midpoints align? sounds like an easyish check
apply(largefish_cumulprop>.5, 2, which.max)
plot(apply(sonar_cumulprop[,-1]>.5, 2, which.max), pch=16, col=4)
points(apply(largefish_cumulprop[,-1]>.5, 2, which.max, pch=16, col=2))

sonar_midpts <- apply(sonar_cumulprop[,-(1:40)]>.5, 2, which.max)
large_midpts <- apply(largefish_cumulprop[,-1]>.5, 2, which.max)

plot(sonar_midpts, ylim=range(sonar_midpts,large_midpts), pch=16, col=4)
points(large_midpts, pch=16, col=2)
plot(sonar_midpts - large_midpts)



# trying a hierarchical logistic regression for PROP LARGE FISH
library(jagshelper)
library(jagsUI)


# bundle data to pass into JAGS
proplarge_data <- list(
  large = round(as.matrix(largefish0[,-1])),
  all = as.matrix(sonar0[,-(1:(ncol(sonar)-ncol(largefish)+1))]),
  day = 1:nrow(largefish),
  nday = nrow(largefish),
  nyear = ncol(largefish)-1
)

# specify model, which is written to a temporary file
proplarge_jags <- tempfile()
cat('model {
  for(j in 1:nyear) {
      large[1,j] ~ dbin(p[1,j], all[1,j])
      trend[1,j] <- b0[j] + b1[j]*day[1]
      resid0[j] ~ dnorm(0, 0.01)
      mu[1,j] <- trend[1,j] + resid0[j]
      logit(p[1,j]) <- mu[1,j]
    for(i in 2:nday) {
      large[i,j] ~ dbin(p[i,j], all[i,j])
      trend[i,j] <- b0[j] + b1[j]*day[i]
      resid[i-1,j] <- mu[i-1,j] - trend[i-1,j]
      # mu[i,j] <- trend[i,j] + phi*resid[i-1,j]
      mu[i,j] ~ dnorm(trend[i,j] + phi*resid[i-1,j], tau)
      logit(p[i,j]) <- mu[i,j]
    }
  }

  for(j in 1:nyear) {
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b1[j] ~ dnorm(mu_b1, tau_b1)
  }

  mu_b0 ~ dnorm(0, 0.001)
  tau_b0 <- pow(sig_b0, -2)
  sig_b0 ~ dunif(0, 10)

  mu_b1 ~ dnorm(0, 0.001)
  tau_b1 <- pow(sig_b1, -2)
  sig_b1 ~ dunif(0, 10)

  phi ~ dunif(0,.9) # dbeta(5,5)
  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)


}', file=proplarge_jags)


# JAGS controls
niter <- 100000
ncores <- 3
# ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  proplarge_jags_out <- jagsUI::jags(model.file=proplarge_jags, data=proplarge_data,
                                     parameters.to.save=c("b0","b1","sig_b0","sig_b1","p","phi","mu","sig"),
                                     n.chains=ncores, parallel=T, n.iter=niter,
                                     n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

# nbyname(proplarge_jags_out)
plotRhats(proplarge_jags_out)
traceworstRhat(proplarge_jags_out, parmfrow = c(3, 3))

par(mfrow=c(2,2))
for(j in 1:proplarge_data$nyear) {
  plot(x=proplarge_data$day,
       y=proplarge_data$large[,j] / proplarge_data$all[,j],
       ylim=c(0, max(proplarge_data$large / proplarge_data$all, na.rm=TRUE)),
       xlab="day", ylab="prop large", type="b")
  envelope(proplarge_jags_out, "p", column=j, add=TRUE)
  curve(expit(proplarge_jags_out$q50$b0[j] +
                proplarge_jags_out$q50$b1[j]*x),
        add=TRUE, lty=2)
  envelope(proplarge_jags_out, "mu", column=j)
  curve(proplarge_jags_out$q50$b0[j] +
                proplarge_jags_out$q50$b1[j]*x,
        add=TRUE, lty=2)
  points(x=proplarge_data$day,
       y=logit(proplarge_data$large[,j] / proplarge_data$all[,j]))
}

par(mfrow=c(2,2))
for(j in 1:proplarge_data$nyear) {
  plot(x=proplarge_data$day,
         y=logit(proplarge_data$large[,j] / proplarge_data$all[,j]),
       type="b")
  abline(lm(logit(proplarge_data$large[,j] / proplarge_data$all[,j]) ~
              proplarge_data$day))
}
