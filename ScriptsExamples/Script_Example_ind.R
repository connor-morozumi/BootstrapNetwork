# load packages, sources and files

require(SYNCSA)
require(picante)
require(bipartite)

setwd("..")
source("./Functions/boot.net.ind.r")

web.ind<-read.table("./DataSetExamples/SS_100individuals.txt", header = TRUE)
web.ind
web <- web.ind[, -1 ]
web
spp <- web.ind[, 1, drop = FALSE]
spp

# To calculate the observed metrics (e.g., NODF) with the interaction matrix. The matrix with the information of the individuals (e.g., "SS_100individuals") is not used here. The observed metric value used in the graphic is from the interaction matrix with the interaction events (e.g., "SS").  
 
webSS <- read.table("./DataSetExamples/SS.txt", header = TRUE) 
webSS

Obs_nodf <- nestednodf(webSS, order = TRUE, weighted = FALSE)$statistic[3]
Obs_nodf

# To calculate the network metric under assessment for the bootstrap samples with replacement method (e.g., NODF)

Res_SS.ind_nodf <- boot.net.ind(web, spp, n.min = 10, by = 5, runs = 1000, method = "nodf", progressbar = TRUE)
Res_SS.ind_nodf

# To calculate the median and confidence limits (lower and upper confidence intervals) of the bootstrap samples 

SS.ind_boot_median <- apply(Res_SS.ind_nodf, 2, median)
SS.ind_boot_median
SS.ind_boot_quantile_lower <- apply(Res_SS.ind_nodf, 2, quantile, probs = 0.025)
SS.ind_boot_quantile_lower
SS.ind_boot_quantile_upper <- apply(Res_SS.ind_nodf, 2, quantile, probs = 0.975)
SS.ind_boot_quantile_upper

# To extract the sample size for the plot

sample.seq <- as.numeric(substr(colnames(Res_SS.ind_nodf), 13, 100)) 
sample.seq

# Plot

plot(SS.ind_boot_median, type = "l", xaxt = "n", ylim = c(0, 100), ylab = "NODF", xlab = "Number of individuals", las = 1) # Draw median
points(SS.ind_boot_quantile_lower, type = "l") # Draw lower quantile
points(SS.ind_boot_quantile_upper, type = "l") # Draw upper quantile
axis(side = 1, at = 1:length(sample.seq), label = c(sample.seq)) # Add axis values and labels
points(length(sample.seq), Obs_nodf, pch = "*", cex = 3) # Add point for the observed nodf
