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


### ---------------- Investigating LARGE FISH PROPORTION ----------------- ###

# calculating large fish proportion
largefish_prop <- largefish #initializing
for(j in 2:ncol(largefish)) {
  largefish_prop[,j] <- largefish[,j]/sonar[[names(largefish)[j]]]
}





# let's see what happens when we remove zeroes
largefishNA <- largefish
largefishNA[largefishNA <= 0] <- NA

# bundle data to pass into JAGS
proplarge_data1 <- list(
  large = round(as.matrix(largefish0[,-1])),
  # large = round(as.matrix(largefishNA[,-1])),
  all = as.matrix(sonar0[,-(1:(ncol(sonar)-ncol(largefish)+1))]),
  day = 1:nrow(largefish),
  day_c = 1:nrow(largefish) - mean(1:nrow(largefish)),  # recentering
  nday = nrow(largefish),
  nyear = ncol(largefish)-1
)
proplarge_data2 <- list(
  # large = round(as.matrix(largefish0[,-1])),
  large = round(as.matrix(largefishNA[,-1])),
  all = as.matrix(sonar0[,-(1:(ncol(sonar)-ncol(largefish)+1))]),
  day = 1:nrow(largefish),
  day_c = 1:nrow(largefish) - mean(1:nrow(largefish)),  # recentering
  nday = nrow(largefish),
  nyear = ncol(largefish)-1
)

# specify model, which is written to a temporary file
model1 <- tempfile()
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
      trend[i,j] <- b0[j] + b1[j]*day_c[i]   # logistic trend wrt day
      # trend[i,j] <- b0[j] + b1[j]*day[i]   # logistic trend wrt day
      mu[i,j] ~ dnorm(trend[i,j], tau)     # tau gives additional overdispersion sd
      logit(p[i,j]) <- mu[i,j]             # logit link
      large[i,j] ~ dbin(p[i,j], all[i,j])  # large fish assumed binomial

      mu_pp[(nday*(j-1)) + i] ~ dnorm(trend[i,j], tau)
    }
  }

  # posterior predictive for a new year with no binomial data
  for(i in 1:nday) {
    trendnew[i] <- b0new + b1new*day_c[i]
    # trendnew[i] <- b0new + b1new*day[i]
    munew[i] ~ dnorm(trendnew[i], tau)
    logit(pnew[i]) <- munew[i]
  }

  for(j in 1:nyear) {
    b0[j] ~ dnorm(mu_b0, tau_b0)
    b1[j] ~ dnorm(mu_b1, tau_b1)
  }
  b0new ~ dnorm(mu_b0, tau_b0)
  b1new ~ dnorm(mu_b1, tau_b1)

  # for(j in 1:nyear) {
  #   b[j, 1:2] ~ dmnorm(mu_b, tau_b)
  #   b0[j] <- b[j,1]
  #   b1[j] <- b[j,2]
  # }
  # bnew  ~ dmnorm(mu_b, tau_b)
  # b0new <- bnew[1]
  # b1new <- bnew[2]
  #
  # tau_b <- inverse(sig_b)
  # sig_b[1,1] <- sig_b0^2
  # sig_b[1,2] <- sig_b0*sig_b1*rho_b
  # sig_b[2,1] <- sig_b0*sig_b1*rho_b
  # sig_b[2,2] <- sig_b1^2
  #
  # mu_b[1] <- mu_b0
  # mu_b[2] <- mu_b1
  #
  # rho_b ~ dunif(-1,1)

  mu_b0 ~ dnorm(0, 0.001)
  tau_b0 <- pow(sig_b0, -2)
  sig_b0 ~ dunif(0, 10)

  mu_b1 ~ dnorm(0, 0.001)
  tau_b1 <- pow(sig_b1, -2)
  sig_b1 ~ dunif(0, 10)

  # phi ~ dunif(0,.9) # dbeta(5,5)
  tau <- pow(sig, -2)
  sig ~ dunif(0, 10)


}', file=model1)

# specify model, which is written to a temporary file
model2 <- tempfile()
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

      mu_pp[(nday*(j-1)) + i] ~ dnorm(trend[i,j], tau)
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


}', file=model2)


# JAGS controls
niter <- 2*1000
# niter <- 100*1000    # 100k in about 2 min on office machine
# ncores <- 3
ncores <- min(10, parallel::detectCores()-1)

outs <- list()
{
  tstart <- Sys.time()
  print(tstart)
  outs[[1]] <- jagsUI::jags(model.file=model1, data=proplarge_data1,
                            parameters.to.save=c("b0","b1","sig_b0","sig_b1","rho_b",
                                                 "p","trend","mu","sig","pnew","trendnew","munew", "mu_pp"),
                            n.chains=ncores, parallel=T, n.iter=niter,
                            n.burnin=niter/2, n.thin=niter/2000)
  outs[[2]] <- jagsUI::jags(model.file=model1, data=proplarge_data2,
                            parameters.to.save=c("b0","b1","sig_b0","sig_b1","rho_b",
                                                 "p","trend","mu","sig","pnew","trendnew","munew", "mu_pp"),
                            n.chains=ncores, parallel=T, n.iter=niter,
                            n.burnin=niter/2, n.thin=niter/2000)
  outs[[3]] <- jagsUI::jags(model.file=model2, data=proplarge_data1,
                            parameters.to.save=c("b0","b1","sig_b0","sig_b1","rho_b",
                                                 "p","trend","mu","sig","pnew","trendnew","munew", "mu_pp"),
                            n.chains=ncores, parallel=T, n.iter=niter,
                            n.burnin=niter/2, n.thin=niter/2000)
  outs[[4]] <- jagsUI::jags(model.file=model2, data=proplarge_data2,
                            parameters.to.save=c("b0","b1","sig_b0","sig_b1","rho_b",
                                                 "p","trend","mu","sig","pnew","trendnew","munew", "mu_pp"),
                            n.chains=ncores, parallel=T, n.iter=niter,
                            n.burnin=niter/2, n.thin=niter/2000)
  print(Sys.time() - tstart)
}

# par(mfrow=c(1,1))
par(mfrow=c(2,2))
for(imodel in 1:4) {
  plotRhats(outs[[imodel]])
}

par(mfrow=c(1,1))
for(imodel in 1:4) {
qq_postpred(outs[[imodel]], p="mu_pp",
            y = logit(as.vector(proplarge_data$large/proplarge_data$all)))
}

for(imodel in 1:4) {
  proplarge_jags_out <- outs[[imodel]]

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
  legend("topright", fill=adjustcolor(c(2,4), 0.5), col=c(2,4), legend=c("trend","pnew"), bty="o")
  envelope(proplarge_jags_out, "trendnew", col=2, add=TRUE, transform="expit")
  envelope(proplarge_jags_out, "munew", main="new year - logit scale", xlab="day")
  envelope(proplarge_jags_out, "trendnew", col=2, add=TRUE)
  legend("topright", fill=adjustcolor(c(2,4), 0.5), col=c(2,4), legend=c("trend","mu"), bty="o")

}


for(imodel in 1:4) {
  proplarge_jags_out <- outs[[imodel]]

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
  # envelope(prop_apportion_allyrs_mcmc,
  #          xlab="Day", ylab="Proportion Large",
  #          main="Logistic model applied to 1980-2018 sonar counts")
  #
  # # adding the years where we actually have prop data
  # for(j in 2:ncol(largefish_cumulprop)) lines(largefish_cumulprop[,j], col=adjustcolor(cols[j],alpha.f=.4), lwd=2)
  # legend("topleft", legend=c("1980-2018 modeled","2019-2025 data"),
  #        fill=c(adjustcolor(4, alpha.f=.5), NA), lwd=c(NA, 2), border=c(4, NA), col=c(NA, "grey"))


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

}

dim(prop_apportion_mcmc_appended)
quantile(prop_apportion_mcmc_appended[,20,], 0.975)
apply(prop_apportion_mcmc_appended[,20,], 2, quantile, 0.975)
apply(prop_apportion_mcmc_appended[,20,], 2, quantile, 0.975) %>%
  median

quantile(prop_apportion_mcmc_appended[,8,], 0.975)
apply(prop_apportion_mcmc_appended[,8,], 2, quantile, 0.975)
apply(prop_apportion_mcmc_appended[,8,], 2, quantile, 0.975) %>%
  median
