## Packages

require(SYNCSA) # for function progressbar
require(picante) # for function matrix2sample
require(bipartite)

## Description

# The boot_net function assesses sampling sufficiency for network metrics, based on bootstrap methods with replacement. It generates confidence intervals for each network metric with increasing sample size. Here we used as examples connectance, NODF and modularity metrics.

## Usage

# boot.net(web, n.min = 30, by = 10, runs = 1000, method = c("nodf", "connectance", "modularity"), progressbar = TRUE)

## Arguments

# web = interaction matrix 
# n.min = minimum number of interaction for bootstrap method with replacement
# by = number of the sample size that will be increased in each step of the bootstrap method
# runs = number of bootstrap samples
# method = "nodf","connectance" or "modularity"
# progressbar = logic argument (TRUE OR FALSE) to show or not the progress bar of the function. Not compatible with RStudio

## Value

# a dataframe with the metrics values, where each column is a sample size (n.min, by, and number maximum of interactions) and row is the result for each bootstrap sample (runs). A bootstrap sample with dimension less than two species in each trophic level will be ingored in the bootstrap step.

boot.net <- function(web, n.min = 30, by = 10, runs = 1000, method = c("nodf", "connectance", "modularity"), progressbar = TRUE){
	if(length(method) > 1){
		stop("\n Only one argument is accepted in method \n")
	}
	if(!(method == "nodf" | method == "connectance" | method == "modularity")){
		stop("\n Invalid method \n")	
	}
	n.interaction <- sum(web)
	if(n.min > n.interaction){
		stop("n.min greater than number of interactions")
	}
	if(!by%%1 == 0){
		stop("by must be an integer")
	}
	if(by > n.interaction){
		stop("by greater than number of interactions")
	}
	n.link <- sum(ifelse(web>0, 1, 0))
	n.row <- dim(web)[1]
	n.col <- dim(web)[2]
	web.1 <- as.matrix(matrix2sample(web))
	web.2 <- matrix(NA, n.interaction, 2)
	k <- 0
	for(i in 1:n.link){
		for(j in 1:as.numeric(web.1[i, 2])){
			k <- k+1
			web.2[k, ] <- c(web.1[i, c(1,3)])
		}
	}
	sample.seq <- seq(n.min, n.interaction, by)
	if(!length(which(sample.seq==n.interaction))>0){
		sample.seq <- c(sample.seq,n.interaction)
	}
	RES <- matrix(NA, runs, length(sample.seq))
	colnames(RES) <- paste("sample.size.", sample.seq, sep = "")
	k <- 0
	l <- 0
	nt <- length(sample.seq)*runs
	for(n in sample.seq){
		k <- k+1
		i <- 0
		while(i<runs){			
			i <- i+1
			l <- l+1
			web.boot <- matrix(0, n.row, n.col)
			colnames(web.boot) <- colnames(web)
			rownames(web.boot) <- rownames(web)			
			sampled <- web.2[sample(1:n.interaction, n, replace = TRUE), ]
			for(j in 1:n){
				row.boot <- which(rownames(web.boot) == sampled[j, 1])
				col.boot <- which(colnames(web.boot) == sampled[j, 2])
				web.boot[row.boot, col.boot] <- web.boot[row.boot, col.boot]+1
			}
			web.boot <- web.boot[!rowSums(web.boot) == 0, !colSums(web.boot) == 0, drop = FALSE]
			if(dim(web.boot)[1]>1 & dim(web.boot)[2]>1){
				if(method == "nodf"){
					RES[i,k] <- nestednodf(web.boot, order = TRUE, weighted = FALSE)$statistic[3]
				}
				if(method == "connectance"){
					RES[i,k] <- networklevel(web.boot, "connectance")
				}		
				if(method == "modularity"){
					log <- capture.output(mod <- computeModules(web.boot, steps = 1E6))
					if(!is.null(mod)){
						RES[i,k] <- mod@likelihood
					}else{
						i <- i-1
						l <- l-1
					}
				}
			}	
			if(progressbar){
				ProgressBAR(l, nt, style = 3)
			}	
		}
	}	
return(RES)
}
