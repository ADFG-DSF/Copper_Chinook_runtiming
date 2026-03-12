### ---------------- loading packages ----------------- ###

library(tidyverse)
library(jagsUI)
library(jagshelper)


### ---------------- loading data ----------------- ###

# STILL NEED TO DEAL WITH DATES FORMATTED AS "9-May" etc

## maybe the better solution is to pivot_longer and ggplot

sonar <- read.csv("flat_data/Copper River sockeye sonar counts 1980-2025.csv")
largefish <- read.csv("flat_data/LargeFish1.csv")
temperature <- read.csv("flat_data/Temperature0.csv") %>%
  filter(`Useable.`==1)





### ---------------- Investigating sonar (ALL FISH) ----------------- ###

# calculating cumulative and cumulative proportion from counts
sonar_cumul <-
  sonar_cumulprop <-
  sonar0 <-
  sonar # initializing
for(j in 2:ncol(sonar)) {
  sonar0[,j] <- ifelse(is.na(sonar[,j]), 0, sonar[,j])  # adding a record where missing counts are treated as 0
  sonar_cumul[,j] <- cumsum(sonar0[,j])
  sonar_cumulprop[,j] <- sonar_cumul[,j]/sum(sonar[,j], na.rm=TRUE)
}

cols <- rcolors(100)

par(mfrow=c(2,2))
plot(NA, main="All sonar targets", ylab="Daily Count", xlab="Day",
     ylim=range(0, sonar[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(sonar)))
for(j in 2:ncol(sonar)) lines(sonar[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(sonar[,-1], na.rm=TRUE), lwd=3)

plot(NA, main="All sonar targets", ylab="Cumulative Count", xlab="Day",
     ylim=range(0, sonar_cumul[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(sonar_cumul)))
for(j in 2:ncol(sonar_cumul)) lines(sonar_cumul[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(sonar_cumul[,-1], na.rm=TRUE), lwd=3)

plot(NA, main="All sonar targets", ylab="Cumulative Proportion", xlab="Day",
     ylim=range(0, sonar_cumulprop[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(sonar_cumulprop)))
for(j in 2:ncol(sonar_cumulprop)) lines(sonar_cumulprop[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(sonar_cumulprop[,-1], na.rm=TRUE), lwd=3)





### ---------------- Investigating LARGE FISH ----------------- ###

# calculating cumulative and cumulative proportion from counts
largefish_cumul <-
  largefish_cumulprop <-
  largefish0 <-
  largefish # initializing
for(j in 2:ncol(largefish)) {
  # adding a record where missing counts are treated as 0
  largefish0[,j] <- ifelse(is.na(largefish[,j]), 0,
                           ifelse(largefish[,j] < 0, 0,
                                  largefish[,j]))
  largefish_cumul[,j] <- cumsum(largefish0[,j])
  largefish_cumulprop[,j] <- largefish_cumul[,j]/sum(largefish[,j], na.rm=TRUE)
}

cols <- rcolors(100)

par(mfrow=c(2,2))
plot(NA, main="Large Fish", ylab="Daily Count", xlab="Day",
     ylim=range(0, largefish[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(largefish)))
for(j in 2:ncol(largefish)) lines(largefish[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(largefish[,-1], na.rm=TRUE), lwd=3)

plot(NA, main="Large Fish", ylab="Cumulative Count", xlab="Day",
     ylim=range(0, largefish_cumul[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(largefish_cumul)))
for(j in 2:ncol(largefish_cumul)) lines(largefish_cumul[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(largefish_cumul[,-1], na.rm=TRUE), lwd=3)

plot(NA, main="Large Fish", ylab="Cumulative Proportion", xlab="Day",
     ylim=range(0, largefish_cumulprop[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(largefish_cumulprop)))
for(j in 2:ncol(largefish_cumulprop)) lines(largefish_cumulprop[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(largefish_cumulprop[,-1], na.rm=TRUE), lwd=3)





### ---------------- Investigating LARGE FISH PROPORTION ----------------- ###

# calculating large fish proportion
largefish_prop <- largefish #initializing
for(j in 2:ncol(largefish)) {
  largefish_prop[,j] <- largefish[,j]/sonar[[names(largefish)[j]]]
}
plot(NA, main="Large Fish Proportion", ylab="Daily Proportion", xlab="Day",
     ylim=range(0, largefish_prop[,-1], na.rm=TRUE), #log="y",
     xlim=c(1, nrow(largefish_prop)))
for(j in 2:ncol(largefish_prop)) lines(largefish_prop[,j], col=adjustcolor(cols[j],alpha.f=.3), lwd=2)
lines(rowMeans(largefish_prop[,-1], na.rm=FALSE), lwd=3)







# do run midpoints align? sounds like an easyish check
sonar_midpts <- apply(sonar_cumulprop[,-(1:40)]>.5, 2, which.max)
large_midpts <- apply(largefish_cumulprop[,-1]>.5, 2, which.max)

plot(sonar_midpts, ylim=range(sonar_midpts,large_midpts), pch=16, col=4)
points(large_midpts, pch=16, col=2)
plot(sonar_midpts - large_midpts)





### ------------------------------------------------------- ###
# trying a hierarchical logistic regression for PROP LARGE FISH
### ------------------------------------------------------- ###

# let's see what happens when we remove zeroes
largefishNA <- largefish
largefishNA[largefishNA <= 0] <- NA

# bundle data to pass into JAGS
proplarge_data <- list(
  # large = round(as.matrix(largefish0[,-1])),
  large = round(as.matrix(largefishNA[,-1])),
  all = as.matrix(sonar0[,-(1:(ncol(sonar)-ncol(largefish)+1))]),
  day = 1:nrow(largefish),
  day_c = 1:nrow(largefish) - mean(1:nrow(largefish)),  # recentering
  nday = nrow(largefish),
  nyear = ncol(largefish)-1
)

# proplarge_data$large <- cbind(round(as.matrix(largefish0[,-1])), 0)
# proplarge_data$all <- cbind(as.matrix(sonar0[,-(1:(ncol(sonar)-ncol(largefish)+1))]), 0)
# proplarge_data$nyear <- ncol(largefish)


# specify model, which is written to a temporary file
proplarge_jags <- tempfile()
cat('model {
  # for(j in 1:nyear) {
  #     large[1,j] ~ dbin(p[1,j], all[1,j])
  #     trend[1,j] <- b0[j] + b1[j]*day[1]
  #     resid0[j] ~ dnorm(0, 0.01)
  #     mu[1,j] <- trend[1,j] + resid0[j]
  #     logit(p[1,j]) <- mu[1,j]
  #   for(i in 2:nday) {
  #     large[i,j] ~ dbin(p[i,j], all[i,j])
  #     trend[i,j] <- b0[j] + b1[j]*day[i]
  #     resid[i-1,j] <- mu[i-1,j] - trend[i-1,j]
  #     # mu[i,j] <- trend[i,j] + phi*resid[i-1,j]
  #     mu[i,j] ~ dnorm(trend[i,j] + phi*resid[i-1,j], tau)
  #     logit(p[i,j]) <- mu[i,j]
  #   }
  # }
  for(j in 1:nyear) {
    for(i in 1:nday) {
      # trend[i,j] <- b0[j] + b1[j]*day_c[i]   # logistic trend wrt day
      trend[i,j] <- b0[j] + b1[j]*day[i]   # logistic trend wrt day
      mu[i,j] ~ dnorm(trend[i,j], tau)     # tau gives additional overdispersion sd
      logit(p[i,j]) <- mu[i,j]             # logit link
      large[i,j] ~ dbin(p[i,j], all[i,j])  # large fish assumed binomial
    }
  }

  # posterior predictive for a new year with no binomial data
  for(i in 1:nday) {
    # trendnew[i] <- b0new + b1new*day_c[i]
    trendnew[i] <- b0new + b1new*day[i]
    munew[i] ~ dnorm(trendnew[i], tau)
    logit(pnew[i]) <- munew[i]
  }

  # for(j in 1:nyear) {
  #   b0[j] ~ dnorm(mu_b0, tau_b0)
  #   b1[j] ~ dnorm(mu_b1, tau_b1)
  # }
  # b0new ~ dnorm(mu_b0, tau_b0)
  # b1new ~ dnorm(mu_b1, tau_b1)

  for(j in 1:nyear) {
    b[j, 1:2] ~ dmnorm(mu_b, tau_b)
    b0[j] <- b[j,1]
    b1[j] <- b[j,2]
  }
  bnew  ~ dmnorm(mu_b, tau_b)
  b0new <- bnew[1]
  b1new <- bnew[2]

  tau_b <- inverse(sig_b)
  sig_b[1,1] <- sig_b0^2
  sig_b[1,2] <- sig_b0*sig_b1*rho_b
  sig_b[2,1] <- sig_b0*sig_b1*rho_b
  sig_b[2,2] <- sig_b1^2

  mu_b[1] <- mu_b0
  mu_b[2] <- mu_b1

  rho_b ~ dunif(-1,1)

  mu_b0 ~ dnorm(0, 0.001)
  tau_b0 <- pow(sig_b0, -2)
  sig_b0 ~ dunif(0, 10)

  mu_b1 ~ dnorm(0, 0.001)
  tau_b1 <- pow(sig_b1, -2)
  sig_b1 ~ dunif(0, 10)

  # phi ~ dunif(0,.9) # dbeta(5,5)
  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)


}', file=proplarge_jags)


# JAGS controls
niter <- 10*1000    # 100k in about 2 min on office machine
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

{
  tstart <- Sys.time()
  print(tstart)
  proplarge_jags_out <- jagsUI::jags(model.file=proplarge_jags, data=proplarge_data,
                                     parameters.to.save=c("b0","b1","sig_b0","sig_b1","rho_b",
                                                          "p","trend","mu","sig","pnew","trendnew","munew"),
                                     n.chains=ncores, parallel=T, n.iter=niter,
                                     n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

## Some diagnostic plots
par(mfrow=c(1,1))
plotRhats(proplarge_jags_out)
traceworstRhat(proplarge_jags_out, parmfrow = c(2, 2))

par(mfrow=c(2,2))
# crossplot(proplarge_jags_out, p=c("b0", "b1"), col="random")
# crossplot(proplarge_jags_out, p=c("b0", "b1"), drawx = TRUE, col="random")
crossplot(proplarge_jags_out, p=c("b0", "b1"), drawblob = TRUE, col="random")
crossplot(proplarge_jags_out, p=c("sig_b0", "sig_b1"), drawblob = TRUE)
crossplot(proplarge_jags_out, p=c("sig", "sig_b0"), drawblob = TRUE)
crossplot(proplarge_jags_out, p=c("sig", "sig_b1"), drawblob = TRUE)

caterpillar(proplarge_jags_out, "b0")
caterpillar(proplarge_jags_out, "b1")
# envelope(proplarge_jags_out, "trend", column=1)
# for(j in 2:7) envelope(proplarge_jags_out, "trend", column=j, add=TRUE)

envelope(proplarge_jags_out, "munew")
for(j in 1:7) curve((proplarge_jags_out$q50$b0[j] + proplarge_jags_out$q50$b1[j]*x), add=TRUE)

envelope(proplarge_jags_out, "pnew")
# for(j in 1:7) lines(proplarge_jags_out$q50$p[,j])
for(j in 1:7) curve(expit(proplarge_jags_out$q50$b0[j] + proplarge_jags_out$q50$b1[j]*x), add=TRUE)

## Plotting model output & data for all years
par(mfrow=c(2,2))
yearnames <- 2019:2025
for(j in 1:proplarge_data$nyear) {
  plot(x=proplarge_data$day,
       y=proplarge_data$large[,j] / proplarge_data$all[,j],
       ylim=c(0, max(proplarge_data$large / proplarge_data$all, na.rm=TRUE)),
       xlab="day", ylab="prop large", type="b", main=yearnames[j])
  legend("topright", fill=adjustcolor(c(2,4), 0.5), col=c(2,4), legend=c("trend","p"), bty="o")
  envelope(proplarge_jags_out, "p", column=j, add=TRUE)
  envelope(proplarge_jags_out, "trend", column=j, col=2, add=TRUE, transform="expit")
  # curve(expit(proplarge_jags_out$q50$b0[j] +
  #               proplarge_jags_out$q50$b1[j]*x),
  #       add=TRUE, lty=2)
  envelope(proplarge_jags_out, "mu", column=j, main=paste(yearnames[j],"- logit scale"), xlab="day")
  envelope(proplarge_jags_out, "trend", column=j, col=2, add=TRUE)
  legend("topright", fill=adjustcolor(c(2,4), 0.5), col=c(2,4), legend=c("trend","mu"), bty="o")
  # curve(proplarge_jags_out$q50$b0[j] +
  #               proplarge_jags_out$q50$b1[j]*x,
  #       add=TRUE, lty=2)
  points(x=proplarge_data$day,
       y=logit(proplarge_data$large[,j] / proplarge_data$all[,j]))
}

## Plotting results for a NEW YEAR
envelope(proplarge_jags_out, "pnew", main="new year", xlab="day", ylab="prop large")
# legend("topright", fill=adjustcolor(c(2,4), 0.5), col=c(2,4), legend=c("trend","pnew"), bty="o")
# envelope(proplarge_jags_out, "trendnew", col=2, add=TRUE, transform="expit")
envelope(proplarge_jags_out, "munew", main="new year - logit scale", xlab="day")
# envelope(proplarge_jags_out, "trendnew", col=2, add=TRUE)
# legend("topright", fill=adjustcolor(c(2,4), 0.5), col=c(2,4), legend=c("trend","mu"), bty="o")




### -------------------------------------------------------------------------- ###
# Applying the p associated with a new year to sonar counts without apportionment!
### -------------------------------------------------------------------------- ###

sonar_toapportion <- sonar[,2:40]  # MAKE THIS MORE ROBUST
sonar_apportion_mcmc <- array(dim=c(dim(proplarge_jags_out$sims.list$pnew),
                                    ncol(sonar_toapportion)))
for(j in 1:ncol(sonar_toapportion)) {
  sonar_apportion_mcmc[,,j] <- proplarge_jags_out$sims.list$pnew *
    matrix(sonar_toapportion[,j],
           nrow=dim(sonar_apportion_mcmc)[1],
           ncol=dim(sonar_apportion_mcmc)[2],
           byrow=TRUE)
}

prop_apportion_mcmc <- cum_apportion_mcmc <- sonar_apportion_mcmc0 <- sonar_apportion_mcmc # initializing
sonar_apportion_mcmc0[is.na(sonar_apportion_mcmc0)] <- 0
for(i in 1:dim(sonar_apportion_mcmc)[1]) {
  for(k in 1:dim(sonar_apportion_mcmc)[3]) {
    cum_apportion_mcmc[i,,k] <- cumsum(sonar_apportion_mcmc0[i,,k])
    prop_apportion_mcmc[i,,k] <- cum_apportion_mcmc[i,,k]/sum(sonar_apportion_mcmc0[i,,k])
  }
}

# combining all years' apportioned prop large
# prop_apportion_allyrs_mcmc <- apply(prop_apportion_mcmc, 1:2, median)  # time-consuming
prop_apportion_allyrs_mcmc <- prop_apportion_mcmc %>%
  apply(., 3, \(x) x, simplify=FALSE) %>%
  do.call(rbind, .)

par(mfrow=c(1,1))
envelope(prop_apportion_allyrs_mcmc,
         xlab="Day", ylab="Proportion Large",
         main="Logistic model applied to 1980-2018 sonar counts")

# adding the years where we actually have prop data
for(j in 2:ncol(largefish_cumulprop)) lines(largefish_cumulprop[,j], col=adjustcolor(cols[j],alpha.f=.4), lwd=2)
legend("topleft", legend=c("1980-2018 modeled","2019-2025 data"),
       fill=c(adjustcolor(4, alpha.f=.5), NA), lwd=c(NA, 2), border=c(4, NA), col=c(NA, "grey"))


# creating and populating an appended proportion array
prop_apportion_mcmc_appended <- array(dim = dim(prop_apportion_mcmc) +
                                        c(0, 0, ncol(largefish)-1))
for(j in 1:dim(prop_apportion_mcmc)[2]) {
  for(k in 1:dim(prop_apportion_mcmc)[3]) {
    prop_apportion_mcmc_appended[, j, k] <- prop_apportion_mcmc[, j, k]
  }
}
for(j in 1:dim(prop_apportion_mcmc)[2]) {
  for(k in (1+dim(prop_apportion_mcmc)[3]):dim(prop_apportion_mcmc_appended)[3]) {
    prop_apportion_mcmc_appended[, j, k] <- largefish_cumulprop[j, k-dim(prop_apportion_mcmc)[3]+1]
  }
}

# combining all years' apportioned prop large
# prop_apportion_allyrs_mcmc_appended <- apply(prop_apportion_mcmc_appended, 1:2, median)
prop_apportion_allyrs_mcmc_appended <- prop_apportion_mcmc_appended %>%
  apply(., 3, \(x) x, simplify=FALSE) %>%
  do.call(rbind, .)

envelope(prop_apportion_allyrs_mcmc_appended,
         xlab="Day", ylab="Proportion Large",
         main="All years aggregated")

# adding the years where we actually have prop data
for(j in 2:ncol(largefish_cumulprop)) lines(largefish_cumulprop[,j], col=adjustcolor(cols[j],alpha.f=.4), lwd=2)
legend("topleft", legend=c("1980-2025 aggregated","2019-2025 data"),
       fill=c(adjustcolor(4, alpha.f=.5), NA), lwd=c(NA, 2), border=c(4, NA), col=c(NA, "grey"))

# overlaying median lines (gets ugly)
med_byyear <- apply(prop_apportion_mcmc_appended, 2:3, median)
for(j in 1:ncol(med_byyear)) {
  lines(med_byyear[,j])
}

# are median envelope bound equivalent to aggregated envelope??
# ANSWER: NO, oddly, they're roughly equivalent to that median across years thing?
u95all <- apply(prop_apportion_mcmc_appended, 2:3, quantile, p=0.975)
l95all <- apply(prop_apportion_mcmc_appended, 2:3, quantile, p=0.025)
u50all <- apply(prop_apportion_mcmc_appended, 2:3, quantile, p=0.75)
l50all <- apply(prop_apportion_mcmc_appended, 2:3, quantile, p=0.25)

u95 <- apply(u95all, 1, median)
l95 <- apply(l95all, 1, median)
u50 <- apply(u50all, 1, median)
l50 <- apply(l50all, 1, median)

# do it the old way
# prop_apportion_allyrs_mcmc_appended <- apply(prop_apportion_mcmc_appended, 1:2, median)

# plot it fresh
envelope(prop_apportion_allyrs_mcmc_appended,
         xlab="Day", ylab="Proportion Large",
         main="All years aggregated")
lines(u95)
lines(l95)
lines(u50)
lines(l50)

# do it the old way
prop_apportion_allyrs_mcmc_appended <- apply(prop_apportion_mcmc_appended, 1:2, median)
