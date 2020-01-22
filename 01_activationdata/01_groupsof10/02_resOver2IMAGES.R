# tables for the paper
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
from_zmap_to_voxelID <- function(nifti.file){
	require(oro.nifti)
	tmp <- readNIfTI(nifti.file)
	out <- which(tmp!=0, arr.ind=TRUE)
	out <-data.frame(out)
	names(out) <- c("X", "Y", "Z")
	out
}


fold <- c("01_N10_smooth_EXCL/")#, "01_N10_smooth6_EXCL/", "02_N20_smooth_EXCL/", "02_N20_smooth6_EXCL/")
p<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/01_smallN/01_results/"
s.p <- "./Onderzoek/doctoraat/reports_abstracts/Paper_finaal/03_manuscript/00_pix/"
for (f in fold){
	# prepare the objects to save the data
	path1 <- paste(p,f, sep="")
	fshort <- substr(f, 4, nchar(f)-1)
	if(which(f == fold) < 3){
		maps <- 500:699	
	}else{
		maps <- 200:399	
	}
	ZZ=maps[1]
	path2 <- paste(path1, sprintf("%04d", ZZ), sep="")
	tmp   <- read.table(paste(path2, "STABOUTVAR.txt", sep="/"), header=TRUE, sep="\t")
	results <- matrix(NA, ncol=10, nrow=1)
	names(results)<- c(names(tmp), "set", "cp")
	for (i in maps){
		path2 <- paste(path1, sprintf("%04d", i), sep="")
		#cat(path2,"\n")
		#tmp <- read.table(paste(path2, "STABOUT75.txt", sep="/"), header=TRUE)
		tmp <- read.table(paste(path2, "STABOUTVAR.txt", sep="/"), header=TRUE)
		tmp2 <- read.table(paste(path2, paste("STABOUTp", sprintf("%04d", i),".txt", sep=""), sep="/"), header=TRUE, sep="\t")
		if(dim(tmp)[1]!= dim(tmp2)[1]) cat(i, " problem \n")
		tmp$set <- i
		tmp$cp <- tmp2$P
		results <- rbind(results, as.matrix(tmp))
	}
	results <- data.frame(results)
	# Determine all Pairwise comparinson of DICE/JACCARD
	OMEGAs <-OVs<-c()
	OMEGAS <-OVS<- c()
	OMEGAfs <-OVfs<-c()
	OMEGAfh <- OVfh<- c()

	n.c.s <- n.c.S <- n.c.fs <-n.c.fh <-  c()
	n.v.s <-n.v.S <- n.v.fs <- n.v.fh <- c()

	### STABILITY SOFT AND HARD
	# all sets -> with  wit voxel id and with voxel coordinates
	cat("stab. Soft and Hard ... \n")
	ALL <- data.frame(clus.id=double(), X=double(), Y=double(), Z=double(), set=double())
	RES.SEL <- rbind(results[which(results$cp<0.2 & results$Stability> 2/3 & results$cp >= 0.01),], results[which(results$cp < 0.01),])
	n.c.s[1] <- mean(table(RES.SEL$set))
	n.c.s[2] <- sd(table(RES.SEL$set))
	n.v.s[1] <- mean(aggregate(size ~ set, data=RES.SEL, sum)[,2])
	n.v.s[2] <- sd(aggregate(size ~ set, data=RES.SEL, sum)[,2])
	SELECTION <- maps
	for (i in SELECTION){
		cat(i, "\t")
		PATH.inA <- paste(path1, sprintf("%04d", i), sep="")
		nfileA <-  paste(PATH.inA, "/cl.nii.gz", sep="")
		idsA <- RES.SEL$ID[which(RES.SEL$set==i)]
		MIa <- min(idsA)
		MAa <- max(idsA)

		tmpA <- as.data.frame(from_cluster_to_voxels(nfileA, MIa, MAa))
		names(tmpA) <- c("clus.id", "X", "Y", "Z")
		tmpA$set <- i
		ALL <- rbind(ALL, tmpA)
	}
	for (k in 1:100){
		i<-SELECTION[((k-1)*2+1)]
		SELECTION2 <-SELECTION[k*2]
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
	cat("FWE. Hard ... \n")
	RES.SEL <- results[which(results$cp<0.01),]
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
		PATH.inA <- paste(path1, sprintf("%04d", i), sep="")
		nfileA <-  paste(PATH.inA, "/cl.nii.gz", sep="")
		idsA <- RES.SEL$ID[which(RES.SEL$set==i)]
		MIa <- min(idsA)
		MAa <- max(idsA)

		tmpA <- as.data.frame(from_cluster_to_voxels(nfileA, MIa, MAa))
		names(tmpA) <- c("clus.id", "X", "Y", "Z")
		tmpA$set <- i
		ALL <- rbind(ALL, tmpA)
	}
	for (k in 1:100){
		i<-SELECTION[((k-1)*2+1)]
		SELECTION2 <-SELECTION[k*2]
		if (length(SELECTION2)==0) next
		A <- ALL[which(ALL$set==i),2:4]
		nA <- dim(A)[1]
	
		omega <-c()
		AB <-c()
		for (j in SELECTION2){
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
	cat("fwe. Soft ... \n")
	RES.SEL <- results[which(results$cp<0.2 ),]
	SELECTION <- unique(RES.SEL$set)
	n.c.fs[1] <- mean(table(RES.SEL$set))
	n.c.fs[2] <- sd(table(RES.SEL$set))
	n.v.fs[1] <- mean(aggregate(size ~ set, data=RES.SEL, sum)[,2])
	n.v.fs[2] <- sd(aggregate(size ~ set, data=RES.SEL, sum)[,2])

	# Determine all Pairwise comparinson of DICE/JACCARD
	# all sets -> with  wit voxel id and with voxel coordinates
	ALL <- data.frame(clus.id=double(), X=double(), Y=double(), Z=double(), set=double())
	for (i in SELECTION){
		#cat(i, "\n")
		PATH.inA <- paste(path1, sprintf("%04d", i), sep="")
		nfileA <-  paste(PATH.inA, "/cl.nii.gz", sep="")
		idsA <- RES.SEL$ID[which(RES.SEL$set==i)]
		MIa <- min(idsA)
		MAa <- max(idsA)

		tmpA <- as.data.frame(from_cluster_to_voxels(nfileA, MIa, MAa))
		names(tmpA) <- c("clus.id", "X", "Y", "Z")
		tmpA$set <- i
		ALL <- rbind(ALL, tmpA)
	}
	for (k in 1:100){
		i<-SELECTION[((k-1)*2+1)]
		SELECTION2 <-SELECTION[k*2]
		if (length(SELECTION2)==0) next
		A <- ALL[which(ALL$set==i),2:4]
		nA <- dim(A)[1]
	
		omega <-c()
		AB <-c()
		for (j in SELECTION2){
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
	cat("stab. Soft only ... \n")
	ALL <- data.frame(clus.id=double(), X=double(), Y=double(), Z=double(), set=double())
	RES.SEL <- results[which(results$cp<0.2 & results$Stability> 2/3),]
	SELECTION <- unique(RES.SEL$set)
	n.c.S[1] <- mean(table(RES.SEL$set))
	n.c.S[2] <- sd(table(RES.SEL$set))
	n.v.S[1] <- mean(aggregate(size ~ set, data=RES.SEL, sum)[,2])
	n.v.S[2] <- sd(aggregate(size ~ set, data=RES.SEL, sum)[,2])
	for (i in SELECTION){
		#cat(i, "\n")
		PATH.inA <- paste(path1, sprintf("%04d", i), sep="")
		nfileA <-  paste(PATH.inA, "/cl.nii.gz", sep="")
		idsA <- RES.SEL$ID[which(RES.SEL$set==i)]
		MIa <- min(idsA)
		MAa <- max(idsA)

		tmpA <- as.data.frame(from_cluster_to_voxels(nfileA, MIa, MAa))
		names(tmpA) <- c("clus.id", "X", "Y", "Z")
		tmpA$set <- i
		ALL <- rbind(ALL, tmpA)
	}
	for (k in 1:100){
		i<-SELECTION[((k-1)*2+1)]
		SELECTION2 <-SELECTION[k*2]
		if (length(SELECTION2)==0) next
		A <- ALL[which(ALL$set==i),2:4]
		nA <- dim(A)[1]
	
		omega <-c()
		AB <-c()
		for (j in SELECTION2){
			#determine the number of voxels in Bj
			B <- ALL[which(ALL$set==j),2:4]
			nB <- dim(B)[1]
			
			#determine the number of voxels in AB
			AuB <-sum(duplicated(rbind(A,B)))
			AB[j] <- AuB
			omega[j] <- AuB/(nA+nB-AuB)
		}
		#lazyness with the indexes from omegad
		OMEGAS <- c(OMEGAS, na.omit(omega))
		OVS <- c(OVS, na.omit(AB))
	}


	dplot <- data.frame(OMEGA=c(OMEGAs, OMEGAS, OMEGAfs, OMEGAfh), correction=c(rep("STAB soft and hard", times=length(OMEGAs)), rep("STAB soft", times=length(OMEGAS)),
	rep("FWE 0.2", times=length(OMEGAfs)), rep("FWE 0.01", times=length(OMEGAfh))))

	out1 <- aggregate(OMEGA~correction, data=dplot, mean)
	out2 <- aggregate(OMEGA~correction, data=dplot, median)
	out3 <- aggregate(OMEGA~correction, data=dplot, sd)

	outO <- cbind(out1, out2[,2], out3[,2])
	names(outO) <- c("scenario", "mean", "median", "sd")

	dplot <- data.frame(OV=c(OVs, OVS, OVfs, OVfh), correction=c(rep("STAB soft and hard", times=length(OVs)), rep("STAB soft", times=length(OVS)),
	rep("FWE 0.2", times=length(OVfs)), rep("FWE 0.01", times=length(OVfh))))
	out1 <- aggregate(OV~correction, data=dplot, mean)
	out2 <- aggregate(OV~correction, data=dplot, median)
	out3 <- aggregate(OV~correction, data=dplot, sd)

	nc <- rbind(n.c.fh, n.c.fs, n.c.S, n.c.s)
	nv <- rbind(n.v.fh, n.v.fs, n.v.S, n.v.s)
	out <- cbind(outO, out1, out2[,2], out3[,2], nc, nv)

	names(out) <- c("scenario", "average $\\omega$", "median $\\omega$", "sd $\\omega$", "mean $n$ overlap", "median $n$ overlap", "sd \\#", "n $c$", "sd $c$", "n $v$ tot", "sd $v$")
	print.xtable(xtable(out, digits=3, type='latex'), sanitize.text.function=identity, file=paste(s.p,fshort, ".tex", sep=""))
}


i=500
path2 <- paste(path1, sprintf("%04d", i), sep="")
tmp   <- read.table(paste(path2, "filecon", sep="/"), header=TRUE, sep="\t")
results <- matrix(NA, ncol=10, nrow=1)
names(results)<- names(tmp)

path2 <- paste(path1, sprintf("%04d", i), sep="")
tmp   <- read.table(paste(path2, "STABOUTVAR.txt", sep="/"), header=TRUE)

results <- matrix(NA, ncol=10, nrow=1)
names(results)<- c(names(tmp), "set", "cp")

for (i in c(500:699)){
	path2 <- paste(path1, sprintf("%04d", i), sep="")
	cat(path2,"\n")
	#tmp <- read.table(paste(path2, "STABOUT75.txt", sep="/"), header=TRUE)
	tmp <- read.table(paste(path2, "STABOUTVAR.txt", sep="/"), header=TRUE)
	tmp2 <- read.table(paste(path2, paste("STABOUTp", sprintf("%04d", i),".txt", sep=""), sep="/"), header=TRUE, sep="\t")
	if(dim(tmp)[1]!= dim(tmp2)[1]) cat(i, " problem \n")
	tmp$set <- i
	tmp$cp <- tmp2$P
	results <- rbind(results, as.matrix(tmp))
}
results<- data.frame(results)
library(xtable)


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
RES.SEL <- rbind(results[which(results$cp<0.2 & results$Stability> 2/3 & results$cp >= 0.01),], results[which(results$cp < 0.01),])
n.c.s[1] <- mean(table(RES.SEL$set))
n.c.s[2] <- sd(table(RES.SEL$set))
n.v.s[1] <- mean(aggregate(size ~ set, data=RES.SEL, sum)[,2])
n.v.s[2] <- sd(aggregate(size ~ set, data=RES.SEL, sum)[,2])

SELECTION <- 500:699
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
for (k in 1:100){
	i<-SELECTION[((k-1)*2+1)]
	SELECTION2 <-SELECTION[k*2]
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







dplot <- data.frame(OMEGA=c(OMEGAs, OMEGAS, OMEGAfs, OMEGAfh), correction=c(rep("STAB soft and hard", times=100), rep("STAB soft", times=100),
rep("FWE 0.2", times=100), rep("FWE 0.01", times=100)))
library(ggplot2)
ggplot(dplot, aes(y=OMEGA, x=correction))+ geom_boxplot()
out1 <- aggregate(OMEGA~correction, data=dplot, mean)
out2 <- aggregate(OMEGA~correction, data=dplot, median)
out3 <- aggregate(OMEGA~correction, data=dplot, sd)

nc <- rbind(n.c.fh, n.c.fs, n.c.S, n.c.s)
nv <- rbind(n.v.fh, n.v.fs, n.v.S, n.v.s)

out <- cbind(out1, out2[,2], out3[,2])
names(out) <- c("scenario", "mean", "median", "sd")

xtable(out, digits=3)

dplot <- data.frame(OV=c(OVs, OVS, OVfs, OVfh), correction=c(rep("STAB soft and hard", times=100), rep("STAB soft", times=100),
rep("FWE 0.2", times=100), rep("FWE 0.01", times=100)))
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


