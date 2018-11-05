# load packages, sources and files

require(SYNCSA)
require(picante)
require(bipartite)

setwd("..")
source("./Functions/boot.net.r")

webSS <- read.table("./DataSetExamples/SS.txt", header = TRUE) 
webSS

# To calculate the observed metrics values 

Obs_nodf <- nestednodf(webSS, order = TRUE, weighted = FALSE)$statistic[3]
Obs_nodf
Obs_connectance <- networklevel(webSS, "connectance")
Obs_connectance
Obs_mod <- computeModules(webSS, steps = 1E6)
Obs_mod <- Obs_mod@likelihood
Obs_mod

# To calculate the network metric under assessment for the bootstrap samples with replacement method (e.g., NODF)

Res_SS_nodf <- boot.net(webSS, n.min = 10, by = 5, runs = 1000, method = "nodf", progressbar = TRUE)
Res_SS_nodf

# To calculate the median and confidence limits (lower and upper confidence intervals) of the bootstrap samples 

SS_boot_median <- apply(Res_SS_nodf, 2, median)
SS_boot_median
SS_boot_quantile_lower <- apply(Res_SS_nodf, 2, quantile, probs = 0.025)
SS_boot_quantile_lower
SS_boot_quantile_upper <- apply(Res_SS_nodf, 2, quantile, probs = 0.975)
SS_boot_quantile_upper

# To extract the sample size for the plot

sample.seq <- as.numeric(substr(colnames(Res_SS_nodf), 13, 100)) 
sample.seq

# Plot

plot(SS_boot_median, type = "l", xaxt = "n", ylim = c(0, 100), ylab = "NODF", xlab = "Number of events", las = 1) # Draw median
points(SS_boot_quantile_lower, type = "l") # Draw lower quantile
points(SS_boot_quantile_upper, type = "l") # Draw upper quantile
axis(side = 1, at = 1:length(sample.seq), label = c(sample.seq)) # Add axis values and labels
points(length(sample.seq), Obs_nodf, pch="*", cex = 3) # Add point for the observed nodf
