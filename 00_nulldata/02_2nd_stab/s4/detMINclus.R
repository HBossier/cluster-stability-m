####################
#### TITLE: output -> min cluster extend to obtain an FWE-corrected cluster p-value of 0.01
#### Contents: 	    calculation using FWHM and approximations in the SPM handbook.
#### Source Files: 	
#### Last Modified: 08/01/2016 
####################

##########################################################################################
# Start with clear working space
rm(list=ls())
##########################################################################################
options(warn=-1)
input <- commandArgs(TRUE)
# [1] min p
# [2] path to mask file
# [3] path to smoothess file
# input <- rep(NA, 3)
#input[1] <- 0.01
#input[2] <- "./Onderzoek/Paper04/02_results/02_emotion/01_smallN/1000/output/mask.nii.gz"
#input[3] <- "./Onderzoek/Paper04/02_results/02_emotion/01_smallN/1000/smoothness1"

#print(input)

#code uit the most legendary spmR (*2009;+2010)
resels <- function(field, mask, fwhm){

	xdim <- dim(field)[1]+2
	ydim <- dim(field)[2]+2
	zdim <- dim(field)[3]+2

	mask <- array(0,dim=c(xdim,ydim,zdim))
	mask[2:(xdim-1),2:(ydim-1),2:(zdim-1)] <- field


	Ex <- Ey <- Ez <- Fxy <- Fxz <- Fyz <- C <- 0

	for (i in 1:xdim){
		#cat("=")
		for (j in 1:ydim){
			for (k in 1:zdim){
				if(mask[i,j,k]==1){
				Ex <- ifelse(mask[i+1,j,k]==1,Ex+1,Ex)
				Ey <- ifelse(mask[i,j+1,k]==1,Ey+1,Ey)
				Ez <- ifelse(mask[i,j,k+1]==1,Ez+1,Ez)

				Fxy <- ifelse(mask[i+1,j,k]==1 && mask[i,j+1,k]==1 && mask[i+1,j+1,k]==1,Fxy+1,Fxy)
				Fxz <- ifelse(mask[i+1,j,k]==1 && mask[i,j,k+1]==1 && mask[i+1,j,k+1]==1,Fxz+1,Fxz)
				Fyz <- ifelse(mask[i,j+1,k]==1 && mask[i,j,k+1]==1 && mask[i,j+1,k+1]==1,Fyz+1,Fyz)

				C <- ifelse(mask[i+1,j,k]==1 && mask[i,j+1,k]==1 && mask[i+1,j+1,k]==1 && mask[i,j,k+1]==1 && mask[i+1,j,k+1]==1 && mask[i,j+1,k+1]==1 && mask[i+1,j+1,k+1]==1,C+1,C)
				}
			}
		}
	}

	P <- sum(mask)
	rx <- 1/fwhm[1]
	ry <- 1/fwhm[2]
	rz <- 1/fwhm[3]

	R0 <- P-(Ex+Ey+Ez)+(Fyz+Fxz+Fxy)-C
	R1 <- (Ex-Fxy-Fxz+C)*rx+(Ey-Fxy-Fyz+C)*ry+(Ez-Fxz-Fyz+C)*rz
	R2 <- (Fxy-C)*rx*ry + (Fxz-C)*rx*rz + (Fyz-C)*ry*rz
	R3 <- C*rx*ry*rz

	Res <- c(R0,R1,R2,R3)
	Res
}
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

mask <- oro.nifti::readNIfTI(input[2])
FWHM <- rep(NA,3)
tmp <- read.table(input[3], skip=16, nrows=1)
FWHM[1] <- tmp[,3]
FWHM[2] <- tmp[,7]
FWHM[3] <- tmp[,11]
u=qnorm(0.999)

res <- resels(array(1, dim=c(91,109,91)), mask=mask, fwhm=FWHM)
COR <- sum( res * spm_ECdensity("Z", u))

# implementatie in de code zoals de formuleringen Worsley (2007, spm handbook) [parameter c] en Hayasaka en Nichols (2003) met betrekking tot benadering vooropstellen.
# FWHM -> full width half maximum to apply on a raw field to obtain a similar smoothness as the one obtained in the image.
# D = dimension -> here 3
# u = first level threshold
# c ~ constant as in Worsley (2007)
pro_S <- function(u=qnorm(0.999), size, FWHM=FWHM,D=3, COR){
	c <- (prod(FWHM)*u^(D/2)*(1-pnorm(u)))/ (spm_ECdensity("Z", u)[4]*gamma(D/2+1))
	prob <- exp(-u*(size/c)^(2/D))
	cprob <- prob*COR
	cprob 
}

# correct with -1 so to set the > threshold in a next stage.
out <- min(which( apply(matrix(1:200), 2, pro_S, u=u, FWHM=FWHM, D=3, COR=COR) < as.numeric(as.character(input[1])))) -1

cat(out)


