####################
#### TITLE: output -> min value + number of selected voxels.	
#### Contents: 	    basic calculation
#### Source Files: 	
#### Last Modified: 19/02/2015 (file was born) 	
####################

##########################################################################################
# Start with clear working space
rm(list=ls())
##########################################################################################
options(warn=-1)
if (Sys.info()[1]=='Linux') {
  	# activate input from command line
	input <- commandArgs(TRUE)
  }else{
  	# activate input from command line
	input <- commandArgs(TRUE)
    # fname
}
print(input)
MIN <- as.numeric(as.character(input[3]))
##########################################################################################
OUT <- matrix(NA, nrow=length(input), ncol=2)
for (i in 1){
	vox=min.index=c()
	IN <- read.table(input[i], header=TRUE, sep="\t")
	idx <- IN$Voxels > MIN
	if(length(idx)!=0){
		vox <- sum(IN$Voxels[idx])
		min.index <- ifelse(is.finite(min(IN$Cluster.Index[idx])), min(IN$Cluster.Index[idx]), 0)
	}else{
		vox=min.index=0
	}
	OUT[i,1] <- min.index
	OUT[i,2] <- vox
}
write.table(OUT, file=paste(input[2],"/cltmp",sep=""), row.names=FALSE, col.names=FALSE)