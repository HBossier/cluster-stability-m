fwhm_x <- rep(0,100)
fwhm_y<- rep(0,100)
fwhm_z<- rep(0,100)
fwhm_z<- rep(0,100)

fwhm <- rep(0,100)
for(i in 1:100){
	tmp <- read.table(paste("./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/01_smallN/01_results/01_smResids/", sprintf("%04d",i), "/smoothness", sprintf("%04d",i), sep=""), skip=16, nrows=1)
	fwhm_x[i] <- tmp[,3]
	fwhm_y[i] <- tmp[,7]
	fwhm_z[i] <- tmp[,11]
}

dplot <- data.frame(FWHM=c(fwhm_x, fwhm_y, fwhm_z), DIM=c(rep("X", 100), rep("Y", 100), rep("Z", 100)), bootstrap.sample=c(rep(1:100, times=3)))
library(ggplot2)
ggplot(dplot, aes(x=bootstrap.sample, y=FWHM, color=DIM)) +geom_point() +stat_smooth() +geom_abline(intercept = 2.73051,slope=0, aes(colour="X"), lty=2)+geom_abline(intercept = 3.04185, slope=0, aes(colour="Y"), lty=2)+geom_abline(intercept = 2.6379, slope=0, aes(colour="Z"), lty=2)+
labs(title = "FWHM (vox) boot. sapmle vs. group analysis")

# compute the intrinsic smoothness of the HCP data based on the original residuals. 
fwhm2 <- rep(0,77)
path1='./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/00_rawData/'
files=read.table(paste(path1, "files.txt", sep=""))
tmp <- read.table(paste("./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/01_smallN/01_results/", sprintf("%04d",i), "/smoothness_UN", sprintf("%04d",i), sep=""))
fwhm_x <- rep(0,77)
fwhm_y <- rep(0,77)
fwhm_z <- rep(0,77)
for (i in 1:77){
	tmp <- read.table(paste(path1,i, "/", files[i,], "/","/MNINonLinear/Results/tfMRI_EMOTION_LR/stats/smoothness_UN", sep=""), skip=16, nrows=1)
	fwhm_x[i] <- tmp[,3]
	fwhm_y[i] <- tmp[,7]
	fwhm_z[i] <- tmp[,11]
}

plot(fwhm2, type="b", col="red", main="FWHMx * FWHMy * FWHMz (in voxels)")
lines(fwhm, col="black", type="b")
abline(h=mean(fwhm), col="black", lty=2)
abline(h=mean(fwhm2), col="red", lty=2)

legend(0,2.45, col=c("red", "black"), lty=1, legend=c(expression(paste("based on ", beta)), "based on original images"))
box(lwd=1.25)

# smoothed images:

path1='./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/00_rawData/'
files=read.table(paste(path1, "files.txt", sep=""))
fwhm_x <- rep(0,77)
fwhm_y<- rep(0,77)
fwhm_z<- rep(0,77)
fwhm_z<- rep(0,77)
dlh <- rep(0,77)

for (i in 1:77){
	tmp <- read.table(paste(path1,i, "/", files[i,], "/","MNINonLinear/Results/tfMRI_EMOTION_LR/tfMRI_EMOTION_LR_hp200_s4.feat/stats/smoothness_SM", sep=""), skip=16, nrows=1)
	fwhm_x[i] <- tmp[,3]
	fwhm_y[i] <- tmp[,7]
	fwhm_z[i] <- tmp[,11]
	dlh[i] <- read.table(paste(path1,i, "/", files[i,], "/","MNINonLinear/Results/tfMRI_EMOTION_LR/tfMRI_EMOTION_LR_hp200_s4.feat/stats/smoothness_SM", sep=""), skip=18, nrows=1)[1,2]
	
}

plot(dlh, type="b", col="red", main="DLH factor in FSL")
lines(fwhm, col="black", type="b")
abline(h=mean(fwhm), col="black", lty=2)
abline(h=mean(dlh), col="red", lty=2)

legend(0,2.45, col=c("red", "black"), lty=1, legend=c(expression(paste("based on ", beta)), "based on original images"))
box(lwd=1.25)
