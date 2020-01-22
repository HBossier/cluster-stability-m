library(oro.nifti)
from_cluster_to_voxels <- function(nifti.file, MI, MA){
	require(oro.nifti)
	tmp <- readNIfTI(nifti.file)
	raCl <- MI:MA
	TOT <- length(which(tmp %in% c(MI:MA)))
	OUT <- matrix(NA, ncol=4, nrow=TOT)
	strt <-1
	for(i in seq(MI, MA)){
		strt1<-length(which(tmp==i))+strt -1
		OUT[strt:strt1,1] <- i
		OUT[strt:strt1,2:4] <- which(tmp==i, arr.ind=TRUE	)
		strt <- strt1+1
	}
	rm(tmp)
	OUT
}

RowMatrix <- function(Row, MATRIX){
	if(length(Row) == dim(MATRIX)[2]){
		out <- c()
		for(i in 1:dim(MATRIX)[1]){
			out[i] <- sum(c(Row[1]==MATRIX[i,1],  Row[2]==MATRIX[i,2], Row[3]==MATRIX[i,3]))==3
		}
	}else{
		warning("check dimensions...")
		out <- NA
	}
	sum(out)
}

#
from_zmap_to_voxelID <- function(nifti.file){
	require(oro.nifti)
	tmp <- readNIfTI(nifti.file)
	out <- which(tmp!=0, arr.ind=TRUE)
	out <-data.frame(out)
	names(out) <- c("X", "Y", "Z")
	out
}

path1 <-p<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/01_smallN/01_results/"
#path1 <- p <-"./Onderzoek/Paper04/02_results/02_emotion/01_smallN/00"

i=1
path2 <- paste(path1, sprintf("%04d", i), sep="")
tmp   <- read.table(paste(path2, "filecon", sep="/"), header=TRUE, sep="\t")
results <- matrix(NA, ncol=10, nrow=1)
names(results)<- names(tmp)

for (i in 1:100){
	path2 <- paste(path1, sprintf("%02d", i), sep="")
	cat(path2,"\n")
	tmp   <- read.table(paste(path2, "filecon", sep="/"), header=TRUE, sep="\t")
	tmp$ID <- i
	results <- rbind(results, as.matrix(tmp))
}

#remove NA row
results <- results[-1,]
results <- as.data.frame(results)

#make data-set with all >75 clusters
rS <- results[which(results$Voxels>41),]
rH <- results[which(results$Voxels>67),]

# What is the amount of clusters that is larger than 10
a   <- table(rS$ID)
a   <- rbind(a,table(rH$ID))
dif <- (a[1,]-a[2,])

# 26/30 from the 30 samples detected clusters that are larger than 10 -> none larger than 75.
# first hard criterion is not reached what to do.

i=1
path2 <- paste(path1, sprintf("%04d", i), sep="")
tmp   <- read.table(paste(path2, "STABOUTVAR.txt", sep="/"), header=TRUE)

results <- matrix(NA, ncol=9, nrow=1)
names(results)<- c(names(tmp), "set")

for (i in c(1:100)){
	path2 <- paste(path1, sprintf("%04d", i), sep="")
	cat(path2,"\n")
	#tmp <- read.table(paste(path2, "STABOUT75.txt", sep="/"), header=TRUE)
	tmp <- read.table(paste(path2, "STABOUTVAR.txt", sep="/"), header=TRUE)
	tmp$set <- i
	results <- rbind(results, as.matrix(tmp))
}
results<- data.frame(results)
library(xtable)
table(results$set[which(results$size>40 & results$Stability> 2/3 & results$size<=67)]) 


# Determine all Pairwise comparinson of DICE/JACCARD
OMEGAs <-c()
OVs<-c()
OMEGAS <- c()
OVS<-c()
OMEGAfs <-c()
OVfs<-c()
OMEGAfh <- c()
OVfh<- c()

n.c.s <- c()
n.c.S <- c()
n.c.fs <- c()
n.c.fh <- c()

n.v.s <- c()
n.v.S <- c()
n.v.fs <- c()
n.v.fh <- c()



### STABILITY SOFT AND HARD
# all sets -> with  wit voxel id and with voxel coordinates
ALL <- data.frame(clus.id=double(), X=double(), Y=double(), Z=double(), set=double())
RES.SEL <- rbind(results[which(results$size>40 & results$Stability> 2/3 & results$size<=67),], results[which(results$size>67),])
n.c.s[1] <- mean(table(RES.SEL$set))
n.c.s[2] <- sd(table(RES.SEL$set))
n.v.s[1] <- mean(aggregate(size ~ set, data=RES.SEL, sum)[,2])
n.v.s[2] <- sd(aggregate(size ~ set, data=RES.SEL, sum)[,2])

SELECTION <- unique(RES.SEL$set)
for (i in SELECTION){
	cat(i, "\t")
	PATH.inA <- paste(p, sprintf("%04d", i), sep="")
	nfileA <-  paste(PATH.inA, "/cl.nii.gz", sep="")
	idsA <- RES.SEL$ID[which(RES.SEL$set==i)]
	MIa <- min(idsA)
	MAa <- max(idsA)

	tmpA <- as.data.frame(from_cluster_to_voxels(nfileA, MIa, MAa))
	names(tmpA) <- c("clus.id", "X", "Y", "Z")
	tmpA$set <- i
	ALL <- rbind(ALL, tmpA)
}
for (i in SELECTION){
	cat(i, "\n")
	SELECTION2 <- SELECTION[which(SELECTION>i)]
	if (length(SELECTION2)==0) next
	A <- ALL[which(ALL$set==i),2:4]
	nA <- dim(A)[1]
	
	omega <-c()
	AB <-c()
	for (j in SELECTION2){
		cat(j, "\t")
		#determine the number of voxels in Bj
		B <- ALL[which(ALL$set==j),2:4]
		nB <- dim(B)[1]
			
		#determine the number of voxels in AB
		AuB <-sum(duplicated(rbind(A,B)))
		AB[j] <- AuB
		omega[j] <- AuB/(nA+nB-AuB)
	}
	#lazyness with the indexes from omega
	OMEGAs <- c(OMEGAs, na.omit(omega))
	OVs <- c(OVs,na.omit(AB))
}


### FWE HARD
# classical FWE p=0.01 NO stability
RES.SEL <- results[which(results$size>67),]
SELECTION2 <- unique(RES.SEL$set)
n.c.fh[1] <- mean(table(RES.SEL$set))
n.c.fh[2] <- sd(table(RES.SEL$set))
n.v.fh[1] <- mean(aggregate(size ~ set, data=RES.SEL, sum)[,2])
n.v.fh[2] <- sd(aggregate(size ~ set, data=RES.SEL, sum)[,2])

# Determine all Pairwise comparinson of DICE/JACCARD
# all sets -> with  wit voxel id and with voxel coordinates
ALL <- data.frame(clus.id=double(), X=double(), Y=double(), Z=double(), set=double())
for (i in SELECTION){
	cat(i, "\t")
	PATH.inA <- paste(p, sprintf("%04d", i), sep="")
	nfileA <-  paste(PATH.inA, "/cl.nii.gz", sep="")
	idsA <- RES.SEL$ID[which(RES.SEL$set==i)]
	MIa <- min(idsA)
	MAa <- max(idsA)

	tmpA <- as.data.frame(from_cluster_to_voxels(nfileA, MIa, MAa))
	names(tmpA) <- c("clus.id", "X", "Y", "Z")
	tmpA$set <- i
	ALL <- rbind(ALL, tmpA)
}

#classical 
for (i in SELECTION){
	cat(i, "\n")
	SELECTION2 <- SELECTION[which(SELECTION>i)]
	if (length(SELECTION2)==0) next
	A <- ALL[which(ALL$set==i),2:4]
	nA <- dim(A)[1]
	
	omega <-c()
	AB <-c()
	for (j in SELECTION2){
		cat(j, "\t")
		#determine the number of voxels in Bj
		B <- ALL[which(ALL$set==j),2:4]
		nB <- dim(B)[1]
			
		#determine the number of voxels in AB
		AuB <-sum(duplicated(rbind(A,B)))
		AB[j] <- AuB
		omega[j] <- AuB/(nA+nB-AuB)
	}
		#lazyness with the indexes from omega
	OMEGAfh <- c(OMEGAfh, na.omit(omega))
	OVfh <- c(OVfh, na.omit(AB))
}

### FWE SOFT
# classical FWE p=0.2 NO stability
RES.SEL <- results[which(results$size>41 ),]
SELECTION <- unique(RES.SEL$set)
n.c.fs[1] <- mean(table(RES.SEL$set))
n.c.fs[2] <- sd(table(RES.SEL$set))
n.v.fs[1] <- mean(aggregate(size ~ set, data=RES.SEL, sum)[,2])
n.v.fs[2] <- sd(aggregate(size ~ set, data=RES.SEL, sum)[,2])

# Determine all Pairwise comparinson of DICE/JACCARD
# all sets -> with  wit voxel id and with voxel coordinates
ALL <- data.frame(clus.id=double(), X=double(), Y=double(), Z=double(), set=double())
for (i in SELECTION){
	cat(i, "\n")
	PATH.inA <- paste(p, sprintf("%04d", i), sep="")
	nfileA <-  paste(PATH.inA, "/cl.nii.gz", sep="")
	idsA <- RES.SEL$ID[which(RES.SEL$set==i)]
	MIa <- min(idsA)
	MAa <- max(idsA)

	tmpA <- as.data.frame(from_cluster_to_voxels(nfileA, MIa, MAa))
	names(tmpA) <- c("clus.id", "X", "Y", "Z")
	tmpA$set <- i
	ALL <- rbind(ALL, tmpA)
}
for (i in SELECTION){
	cat(i, "\n")
	SELECTION2 <- SELECTION[which(SELECTION>i)]
	if (length(SELECTION2)==0) next
	A <- ALL[which(ALL$set==i),2:4]
	nA <- dim(A)[1]
	omega <-c()
	AB <-c()
	for (j in SELECTION2){
		cat(j, "\t")
		#determine the number of voxels in Bj
		B <- ALL[which(ALL$set==j),2:4]
		nB <- dim(B)[1]
			
		#determine the number of voxels in AB
		AuB <-sum(duplicated(rbind(A,B)))
		AB[j] <- AuB
		omega[j] <- AuB/(nA+nB-AuB)
	}
		#lazyness with the indexes from omega
	OMEGAfs <- c(OMEGAfs, na.omit(omega))
	OVfs <- c(OVfs, na.omit(AB))
}

### STABILITY SOFT 
# all sets -> with  wit voxel id and with voxel coordinates
ALL <- data.frame(clus.id=double(), X=double(), Y=double(), Z=double(), set=double())
RES.SEL <- results[which(results$size>40 & results$Stability> 2/3),]
SELECTION <- unique(RES.SEL$set)
n.c.S[1] <- mean(table(RES.SEL$set))
n.c.S[2] <- sd(table(RES.SEL$set))
n.v.S[1] <- mean(aggregate(size ~ set, data=RES.SEL, sum)[,2])
n.v.S[2] <- sd(aggregate(size ~ set, data=RES.SEL, sum)[,2])

for (i in SELECTION){
	cat(i, "\n")
	PATH.inA <- paste(p, sprintf("%04d", i), sep="")
	nfileA <-  paste(PATH.inA, "/cl.nii.gz", sep="")
	idsA <- RES.SEL$ID[which(RES.SEL$set==i)]
	MIa <- min(idsA)
	MAa <- max(idsA)

	tmpA <- as.data.frame(from_cluster_to_voxels(nfileA, MIa, MAa))
	names(tmpA) <- c("clus.id", "X", "Y", "Z")
	tmpA$set <- i
	ALL <- rbind(ALL, tmpA)
}
for (i in SELECTION){
	cat(i, "\n")
	SELECTION2 <- SELECTION[which(SELECTION>i)]
	if (length(SELECTION2)==0) next
	A <- ALL[which(ALL$set==i),2:4]
	nA <- dim(A)[1]
	
	omega <-c()
	AB <-c()
	for (j in SELECTION2){
		cat(j, "\t")
		#determine the number of voxels in Bj
		B <- ALL[which(ALL$set==j),2:4]
		nB <- dim(B)[1]
			
		#determine the number of voxels in AB
		AuB <-sum(duplicated(rbind(A,B)))
		AB[j] <- AuB
		omega[j] <- AuB/(nA+nB-AuB)
	}
	#lazyness with the indexes from omega
	OMEGAS <- c(OMEGAS, na.omit(omega))
	OVS <- c(OVS, na.omit(AB))
}




4656
4950
dplot <- data.frame(OMEGA=c(OMEGAs, OMEGAS, OMEGAfs, OMEGAfh), correction=c(rep("STAB soft and hard", times=4950), rep("STAB soft", times=4950),
rep("FWE 0.2", times=4950), rep("FWE 0.01", times=4950)))
library(ggplot2)
ggplot(dplot, aes(y=OMEGA, x=correction))+ geom_boxplot()
out1 <- aggregate(OMEGA~correction, data=dplot, mean)
out2 <- aggregate(OMEGA~correction, data=dplot, median)
out3 <- aggregate(OMEGA~correction, data=dplot, sd)

nc <- rbind(n.c.fh, n.c.fs, n.c.S, n.c.s)
nv <- rbind(n.v.fh, n.v.fs, n.v.S, n.v.s)

out <- cbind(out1, out2[,2], out3[,2])
names(out) <- c("scenario", "mean", "median", "sd")

xtable(out, digits=5)

dplot <- data.frame(OV=c(OVs, OVS, OVfs, OVfh), correction=c(rep("STAB soft and hard", times=4950), rep("STAB soft", times=4950),
rep("FWE 0.2", times=4950), rep("FWE 0.01", times=4950)))
out1 <- aggregate(OV~correction, data=dplot, mean)
out2 <- aggregate(OV~correction, data=dplot, median)
out3 <- aggregate(OV~correction, data=dplot, sd)
summary(lm(OV ~ correction, dplot))

nc <- rbind(n.c.fh, n.c.fs, n.c.S, n.c.s)
nv <- rbind(n.v.fh, n.v.fs, n.v.S, n.v.s)
out <- cbind(out1, out2[,2], out3[,2], nc, nv)

names(out) <- c("scenario", "mean", "median", "sd", "n $c$", "sd $c$", "n $v$", "sd $v$")
xtable(out, digits=3)



save(OMEGA, file="omega.Rdata")


# Compare with voxel-wise FWE 0.05 over ALL subjects in the 80 group.
ref <- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/01_analysis/01_SMOOTH/results/thresholded.nii.gz"
SELECTION_VOX <- from_zmap_to_voxelID(ref)

OMEGA_withVOX<- rep(NA, 50)
for (i in SELECTION){
	cat(i, "\n")
	A <- ALL[which(ALL$set==i),2:4]
	nA <- dim(A)[1]
		
	B <- SELECTION_VOX
	nB <- dim(B)[1]
	#determine the number of voxels in AB
	AuB <-sum(duplicated(rbind(A,B)))
	OMEGA_withVOX[i] <- AuB/(nA+nB-AuB)
}

# Voxelwise comparison

path1 <-p<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/01_smallN/01_results/00"
i=1
path2 <- paste(path1, sprintf("%02d/", i), "thresholded.nii.gz", sep="")

test <- from_zmap_to_voxelID(path2)
