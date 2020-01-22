# Code to add a $p$ -value to the sub-threshold clusters and to determine the over-all FWE rates
# in the resting-state NULL data. 

# code from the most legendary spmR (*2009;+2010) package to determine Euler Characteristic Density.
spm_ECdensity <- function(STAT,t,df){

      t <- t(t)
      EC <- matrix(rep(NA,4),nrow=4)

      if (STAT=="Z"){
      a <- 4*log(2)
      b <- exp(-t^2/2)

      EC[1,] <- 1-pnorm(t)
      EC[2,] <- a^(1/2)/(2*pi)*b
      EC[3,] <- a/((2*pi)^(3/2))*b*t
      EC[4,] <- a^(3/2)/((2*pi)^2)*b*(t^2-1)
      }
      if (STAT=="T"){
      v <- df[2]
      a <- 4*log(2)
      b <- exp(log(gamma((v+1)/2)) - log(gamma(v/2)))
      c <- (1+t^2/v)^((1-v)/2)

      EC[1,] <- 1-pt(t,v)
      EC[2,] <- a^(1/2)/(2*pi)*c
      EC[3,] <- a/((2*pi)^(3/2))*c*t/((v/2)^(1/2))*b
      EC[4,] <- a^(3/2)/((2*pi)^2)*c*((v-1)*(t^2)/v-1)
      }

      if (STAT=="X"){
      v <- df[2]
      a <- (4*log(2))/(2*pi)
      b <- t^(1/2*(v-1))*exp(-t/2-log(gamma(v/2)))/2^((v-2)/2)

      EC[1,] <- 1-pchisq(t,v)
      EC[2,] <- a^(1/2)*b
      EC[3,] <- a*b*(t-(v-1))
      EC[4,] <- a^(3/2)*b*(t^2-(2*v-1)*t+(v-1)*(v-2))
      }

      if (STAT=="F"){
      k <- df[2]
      v <- df[1]
      a <- (f*log(2))/(2*pi)
      b <- log(gamma(v/2)) + log(gamma(k/2))

      EC[1,] <- 1-pf(t,df)
      EC[2,] <- a^(1/2)*exp(log(gamma((v+k-1)/2))-b)*2^(1/2)*(k*t/v)^(1/2*(k-1))*(1+k*t/v)^(-1/2*(v+k-2))
      EC[3,] <- a*exp(log(gamma((v+k-2)/2))-b)*(k*t/v)^(1/2*(k-2))*(1+k*t/v)^(-1/2*(v+k-2))*((v-1)*k*t/v-(k-1))
      EC[4,] <- a^(3/2)*exp(log(gamma(((v+k-3)/2)))-b)*2^(-1/2)*(k*t/v)^(1/2*(k-3))*(1+k*t/v)^(-1/2*(v+k-2))*((v-1)*(v-2)*(k*t/v)^2-2(2*v*k-v-k-1)*(k*t/v)+(k-1)*(k-2))
      }
      EC
}

# implementatie in de code zoals de formuleringen Worsley (2007, spm handbook) [parameter c] en Hayasaka en Nichols (2003) met betrekking tot benadering vooropstellen.
# FWHM -> full width half maximum to apply on a raw field to obtain a similar smoothness as the one obtained in the image.
# D = dimension -> here 3
# u = first level threshold
# c ~ constant as in Worsley (2007)
pro_S <- function(u=qnorm(0.999), size, FWHM=AV_FWHM,D=3){
	c <- (prod(FWHM)*u^(D/2)*(1-pnorm(u)))/ (spm_ECdensity("Z", u)[4]*gamma(D/2+1))
	prob <- exp(-u*(size/c)^(2/D))
	prob
}

p.tmp<-"./Onderzoek/Paper04/01_data/00_null_s6"

FWE001 	<- data.frame(avsize=double(), a001=double(), as0016=double(),as0017=double(), as00167=double(), as0018=double(),as0019=double())
FWE005	<- data.frame(avsize=double(), a005=double(), as0056=double(),as0057=double(), as00567=double(), as0058=double(),as0019=double())
FWE01	<- data.frame(avsize=double(), a01=double(), as016=double(), as017=double(), as0167=double(), as018=double(),as019=double())
FWE05 	<- data.frame(avsize=double(), a05=double(), as056=double(), as057=double(), as0567=double(), as058=double(),as059=double())
for(i in 1:1000){
if(i%%100==0)cat(i/10, "% .. ")
# check wether cluster are detected:
control.file <- paste(p.tmp, "/", sprintf("%04d", i),"/","filecon", sep="")
control.file2 <- paste(p.tmp, "/", sprintf("%04d", i),"/","smoothness",sprintf("%04d", i), sep="")
control.file3 <- paste(p.tmp, "/", sprintf("%04d", i),"/","STAB_OUT001.txt", sep="")
control.file4 <- paste(p.tmp, "/", sprintf("%04d", i),"/","STAB_OUT005.txt", sep="")
control.file5 <- paste(p.tmp, "/", sprintf("%04d", i),"/","STAB_OUT01.txt", sep="")
control.file6 <- paste(p.tmp, "/", sprintf("%04d", i),"/","STAB_OUT05.txt", sep="")

if(file.exists(control.file) & file.exists(control.file2) & file.exists(control.file3) & file.exists(control.file4) & file.exists(control.file5) & file.exists(control.file6)){
tmp <-read.table(control.file, sep="\t", header=TRUE)
if(dim(tmp)[1]==0){
FWE001[i,]<- 0
FWE005[i,]<- 0
FWE01[i,]<- 0
FWE05[i,]<- 0
}
else{
# file paths
sm.path <- paste(p.tmp, "/", sprintf("%04d", i),"/","smoothness",sprintf("%04d", i), sep="")
fc.path001 <- paste(p.tmp, "/", sprintf("%04d", i),"/","STAB_OUT001.txt", sep="")
fc.path005 <- paste(p.tmp, "/", sprintf("%04d", i),"/","STAB_OUT005.txt", sep="")
fc.path01 <- paste(p.tmp, "/", sprintf("%04d", i),"/","STAB_OUT01.txt", sep="")
fc.path05 <- paste(p.tmp, "/", sprintf("%04d", i),"/","STAB_OUT05.txt", sep="")

# get all information from the output from fslsmoothest with -V option on. 
FWHM <- read.table(sm.path, skip=16, nrows=1)[c(3,7,11)]
resels <- read.table(sm.path, skip=19, nrows=1)[c(2)]/read.table(sm.path, skip=20, nrows=1)[c(2)]
# set the first threhsold
u <- qnorm(0.999)
# compute the correction factor
COR <- sum(unlist(resels)*spm_ECdensity("Z", u))

sizes001 <- read.table(fc.path001, header=TRUE)#; cat(dim(sizes001))
sizes005 <- read.table(fc.path005, header=TRUE)
sizes01 <- read.table(fc.path01, header=TRUE)
sizes05 <- read.table(fc.path05, header=TRUE)

sizes001$p <- COR*apply(as.matrix(sizes001$size),c(1,2), FUN=pro_S, u=qnorm(0.999), FWHM=FWHM, D=3)
sizes005$p <- COR*apply(as.matrix(sizes005$size),c(1,2), FUN=pro_S, u=qnorm(0.999), FWHM=FWHM, D=3)
sizes01$p  <- COR*apply(as.matrix(sizes01$size) ,c(1,2), FUN=pro_S, u=qnorm(0.999), FWHM=FWHM, D=3)
sizes05$p  <- COR*apply(as.matrix(sizes05$size) ,c(1,2), FUN=pro_S, u=qnorm(0.999), FWHM=FWHM, D=3)

FWE001[i,]<- c(mean(sizes001[which(sizes001$p<0.001),"size"], na.rm=TRUE),
dim(sizes001[which(sizes001$p<0.001),])[1],
dim(sizes001[which(sizes001$p<0.001),])[1]+dim(sizes001[which(sizes001$p>0.001 & sizes001$p<0.05 & sizes001$Stability>0.6),])[1],
dim(sizes001[which(sizes001$p<0.001),])[1]+dim(sizes001[which(sizes001$p>0.001 & sizes001$p<0.05 & sizes001$Stability>(2/3)),])[1],
dim(sizes001[which(sizes001$p<0.001),])[1]+dim(sizes001[which(sizes001$p>0.001 & sizes001$p<0.05 & sizes001$Stability>0.7),])[1],
dim(sizes001[which(sizes001$p<0.001),])[1]+dim(sizes001[which(sizes001$p>0.001 & sizes001$p<0.05 & sizes001$Stability>0.8),])[1],
dim(sizes001[which(sizes001$p<0.001),])[1]+dim(sizes001[which(sizes001$p>0.001 & sizes001$p<0.05 & sizes001$Stability>0.9),])[1])

FWE005[i,]<- c(mean(sizes005[which(sizes005$p<0.005),"size"], na.rm=TRUE),
dim(sizes005[which(sizes005$p<0.005),])[1],
dim(sizes005[which(sizes005$p<0.005),])[1]+ dim(sizes005[which(sizes005$p>0.005 & sizes005$p<0.1 & sizes005$Stability>0.6),])[1],
dim(sizes005[which(sizes005$p<0.005),])[1]+ dim(sizes005[which(sizes005$p>0.005 & sizes005$p<0.1 & sizes005$Stability>(2/3)),])[1],
dim(sizes005[which(sizes005$p<0.005),])[1]+ dim(sizes005[which(sizes005$p>0.005 & sizes005$p<0.1 & sizes005$Stability>0.7),])[1],
dim(sizes005[which(sizes005$p<0.005),])[1]+ dim(sizes005[which(sizes005$p>0.005 & sizes005$p<0.1 & sizes005$Stability>0.8),])[1],
dim(sizes005[which(sizes005$p<0.005),])[1]+ dim(sizes005[which(sizes005$p>0.005 & sizes005$p<0.1 & sizes005$Stability>0.9),])[1])

FWE01[i,]<- c(mean(sizes01[which(sizes01$p<0.01),"size"], na.rm=TRUE),
dim(sizes01[which(sizes01$p<0.01),])[1],
dim(sizes01[which(sizes01$p<0.01),])[1]+dim(sizes01[which(sizes01$p>0.01 & sizes01$p<0.2 & sizes01$Stability>0.6),])[1],
dim(sizes01[which(sizes01$p<0.01),])[1]+dim(sizes01[which(sizes01$p>0.01 & sizes01$p<0.2 & sizes01$Stability>(2/3)),])[1],
dim(sizes01[which(sizes01$p<0.01),])[1]+dim(sizes01[which(sizes01$p>0.01 & sizes01$p<0.2 & sizes01$Stability>0.7),])[1],
dim(sizes01[which(sizes01$p<0.01),])[1]+dim(sizes01[which(sizes01$p>0.01 & sizes01$p<0.2 & sizes01$Stability>0.8),])[1],
dim(sizes01[which(sizes01$p<0.01),])[1]+dim(sizes01[which(sizes01$p>0.01 & sizes01$p<0.2 & sizes01$Stability>0.9),])[1])

FWE05[i,]<- c(mean(sizes05[which(sizes05$p<0.05),"size"], na.rm=TRUE),
dim(sizes05[which(sizes05$p<0.05),])[1],
dim(sizes05[which(sizes05$p<0.05),])[1]+dim(sizes05[which(sizes05$p>0.05 & sizes05$p<0.4 & sizes05$Stability>0.6),])[1],
dim(sizes05[which(sizes05$p<0.05),])[1]+dim(sizes05[which(sizes05$p>0.05 & sizes05$p<0.4 & sizes05$Stability>(2/3)),])[1],
dim(sizes05[which(sizes05$p<0.05),])[1]+dim(sizes05[which(sizes05$p>0.05 & sizes05$p<0.4 & sizes05$Stability>0.7),])[1],
dim(sizes05[which(sizes05$p<0.05),])[1]+dim(sizes05[which(sizes05$p>0.05 & sizes05$p<0.4 & sizes05$Stability>0.8),])[1],
dim(sizes05[which(sizes05$p<0.05),])[1]+dim(sizes05[which(sizes05$p>0.05 & sizes05$p<0.4 & sizes05$Stability>0.9),])[1])


}
# if something went wrong with the cluster algorithm
}
else{
FWE001[i,]<- NA
FWE005[i,]<- NA
FWE01[i,]<- NA
FWE05[i,]<- NA
}
}
# NO FP with WLS estimation from FLAME1 in FSL
save(FWE001, FWE005, FWE01, FWE05, file="./Onderzoek/Paper04/01_data/03_null/s6_10_FWE.Rdata")
# what about OLS?
