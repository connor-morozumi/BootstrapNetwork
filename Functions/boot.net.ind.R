## Packages

require(SYNCSA) # for function progressbar
require(picante) # for function matrix2sample
require(bipartite)

# The boot.net function assesses sampling sufficiency for network metrics, based on bootstrap methods with replacement. It generates confidence intervals for each network metric with increasing sample size. Here we used as examples connectance, NODF and modularity metrics, and the number of individuals as sample size. 

## Usage

# boot.net.ind(web, spp, n.min=30, by=10, runs=1000, method=c("nodf", "connectance", "modularity"), progressbar=TRUE)

## Arguments

# web = interaction matrix 
# spp = matrix with individuals code of species, with the same number of lines that the column with the names of individuals. The names of rows must be the code of individuals and the first column the code of species names.
# n.min = minimum number of interaction for bootstrap method with replacement
# by = number of the sample size that will be increased in each step of the bootstrap method
# runs = number of times that the step will be repeated 
# method = "nodf","connectance" or "modularity"
# progressbar = logic argument (TRUE OR FALSE) to show or not the progress bar of the function. Not compatible with RStudio

## Value

# a dataframe with the metrics values, where each column is a sample size (n.min, by, and number maximum of interactions) and row is the result for each bootstrap sample (runs).

boot.net.ind <- function(web, spp, n.min = 30, by = 10, runs = 1000, method = c("nodf", "connectance", "modularity"), progressbar = TRUE){
	if(length(method) > 1){
		stop("\n Only one argument is accepted in method \n")
	}
	if(!(method == "nodf" | method == "connectance" | method == "modularity")){
		stop("\n Invalid method \n")	
	}
	names.spp <- unique(spp[,1])
	n.spp <- length(spp[,1])
	if(n.min > n.spp){
		stop("n.min greater than number of species")
	}
	if(!by%%1 == 0){
		stop("by must be an integer")
	}
	if(by > n.spp){
		stop("by greater than number of species")
	}
	n.col <- dim(web)[2]
	sample.seq <- seq(n.min, n.spp, by)
	if(!length(which(sample.seq == n.spp)) > 0){
		sample.seq <- c(sample.seq, n.spp)
	}
	RES <- matrix(NA, runs, length(sample.seq))
	colnames(RES) <- paste("sample.size.", sample.seq, sep = "")
	k <- 0
	l <- 0
	nt <- length(sample.seq)*runs
	for(n in sample.seq){
		k <- k+1
		i <- 0
		while(i < runs){		
			i <- i+1
			l <- l+1
			web.boot <- matrix(0, length(names.spp), n.col)
			colnames(web.boot) <- colnames(web)
			rownames(web.boot) <- names.spp
			sampled <- sample(rownames(web), n, replace = TRUE)
			for(j in 1:n){
				row.boot <- which(rownames(web.boot) == spp[sampled[j], ])
				web.boot[row.boot, 1:n.col] <- as.numeric(web.boot[row.boot, ]+web[sampled[j], ])
			}
			web.boot <- web.boot[!rowSums(web.boot) == 0, !colSums(web.boot) == 0, drop = FALSE]
			if(method == "nodf"){
				RES[i, k] <- nestednodf(web.boot, order = TRUE, weighted = FALSE)$statistic[3]
			}
			if(method == "connectance"){
				RES[i, k] <- networklevel(web.boot, "connectance")
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
			if(progressbar){
				ProgressBAR(l, nt, style = 3)
			}	
		}
	}	
return(RES)
}