# Structuur Matlab file
# For a single axial slice (Z = 2 mm) of data, you should have:
# Bmap_N_S: 'Difference between Novel and Standard betas averaged over 28 subjects'
# Tmap_N_S: 'T-statistics for the paired t-test comparing Novel and Standard betas'
# Pmap_N_S: 'Binary map indicating significance at P<0.001 (fdr corrected)'
# Underlay: 'Structural image ch2bet from MRIcron, warped to functional data'
rm(list=ls())
library(oro.nifti)
library(R.matlab)
s.path  <- "./Onderzoek/doctoraat/scripts/Paper04/01_data/00_singleSubject/02_results/"
s.path2 <- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_ALL/results/"
r.path <-  "./Onderzoek/doctoraat/scripts/Paper04/01_data/00_singleSubject/01_origData/01/out/"

#### LET OP FSL START OP VOXEL ID == 00 !!! C codering en NIET op 1!
#Low RES
B1 <- readNIfTI(paste(r.path, "out1.nii.gz", sep=""))
Z1 <- readNIfTI(paste(r.path, "zstat1LR.nii.gz", sep=""))
S1 <- readNIfTI(paste(s.path, "OUT-con1.nii.gz", sep=""))
A1 <- readNIfTI(paste(s.path, "MNI152_T1_2mm_brain.nii.gz", sep=""))[,,,1]
B2 <- readNIfTI(paste(s.path2,"out80-con1.nii.gz", sep=""))
Z2 <- readNIfTI(paste(s.path2,"zstat1.nii.gz", sep=""))

z <- Z1[,,53]
b <- B1[,,53]
a <- A1[,,53]
s <- S1[,,53]
b2 <- B2[,,53]
z2 <- Z2[,,53]

filename <- paste(s.path, "CON53.mat", sep="")
writeMat(filename, A=a, B=b, S=s, Z=z, B2=b2)

#Low RES
z <- Z1[,,35]
b <- B1[,,35]
a <- A1[,,35]
s <- S1[,,35]
b2 <- B2[,,35]
z2 <- Z2[,,35]

filename <- paste(s.path, "CON34.mat", sep="")
writeMat(filename, A=a, B=b, S=s, Z=z, B2=b2)

#Low RES
z <- Z1[,,34]
b <- B1[,,34]
a <- A1[,,34]
s <- S1[,,34]
b2 <- B2[,,34]
z2 <- Z2[,,34]

filename <- paste(s.path, "CON33.mat", sep="")
writeMat(filename, A=a, B=b, S=s, Z=z, B2=b2)
