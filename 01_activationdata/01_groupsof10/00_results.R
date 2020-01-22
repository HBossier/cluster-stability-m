# functions
from_cluster_to_voxels <- function(nifti.file){
	require(oro.nifti)
	tmp <- readNIfTI(nifti.file)
	raCl <- range(tmp)
	TOT <- length(which(tmp!=0))
	OUT <- matrix(NA, ncol=4, nrow=TOT)
	strt <-1
	for(i in seq(raCl[1]+1, raCl[2])){
		strt1<-length(which(tmp==i))+strt -1
		#cat(strt, " to ")
		#cat(strt1, "\n")
		OUT[strt:strt1,1] <- i
		OUT[strt:strt1,2:4] <- which(tmp==i, arr.ind=TRUE	)
		strt <- strt1+1
	}
	rm(tmp)
	OUT
}
library(Hmisc)
library(oro.nifti)
cat("start analysing results \n")
#from cluster to voxel
fs <- c("OUT75", "OUTVAR")
for(K in fs){
	nfile 		<- "cl.nii.gz"
	stab.file 	<- paste(K, ".nii.gz", sep="")
	fcon		<- "filecon"
	cat(stab.file, " .. read in .. \n")
	voxMatrix <- from_cluster_to_voxels(nfile)
	tmp<-readNIfTI(stab.file)
	voxMatrix <- cbind(voxMatrix,tmp[voxMatrix[,2:4]])

	voxM <- data.frame(voxMatrix)
	names(voxM) <- c("ID", "x", "y","z", "stab")
	voxM$ID <- as.factor(voxM$ID)
	clS		<- aggregate(stab ~ ID, data=voxM, mean)
	clS_sd 	<- aggregate(stab ~ ID, data=voxM, sd)
	names(clS) <- c("ID", "Stability")
	names(clS_sd) <- c("ID", "sd")
	clS$sd <- clS_sd$sd
	tmp <- readNIfTI(nfile)
	tmpseq <- as.numeric(as.character(clS$ID))
	for(i in tmpseq){
		clS$size[tmpseq==i]<-length(which(tmp==i))
	}
	rm(tmp)

	# Read in the orignal cluster file
	tmp <-read.table(fcon, sep="\t", header=TRUE)
	names(tmp)[1] <- "ID"
	out <- merge(x = clS, y = tmp, by ="ID", all.x=TRUE)

	out<- out[,-c(5,10:12)]
	out<-out[order(out$ID, decreasing=TRUE),]
	names(out)[5:8] <- c("max Z", "X", "Y", "Z")
	out$ID <- as.integer(as.character(out$ID))
	cat("\t", "write output", "\n")
	for(i in out$ID[1:5]){
		pdf(file=paste("cluster_", i, "_", K, ".pdf",sep=""))
		hist(voxM$stab[voxM$ID==i], freq=FALSE, main=paste("cluster", i), xlab="re-selection rate", ylab="density")
		dev.off()
	}

	z <- latex(file=paste("t_emotionSEL", K, ".tex", sep=""),out[which(out$size>75),], digits=3, rowname=NULL, 
	caption="Cluster, ID, stability and SD of stability, overlap with group map and with other image of the SELECTED voxels.", label="T1.emotionS")
	out2 <- out[which(out$size<75 & out$size > 10), ]
	out2<-out2[order(out2$Stability, decreasing=TRUE),]

	z <- latex(file=paste("t_emotionNON-SEL", K, ".tex", sep=""),out2, digits=3, rowname=NULL, 
	caption="Cluster, ID, stability and SD of stability, overlap with group map and with other image of the NON-SELECTED voxels", label="T1.emotionNS")
	write.table(out, file=paste("STAB", K,".txt", sep=""), row.names=FALSE)
	cat("\t . \n")
}
cat("done. \n")