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
  largefish0[,j] <- ifelse(is.na(largefish[,j]), 0, largefish[,j])
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
