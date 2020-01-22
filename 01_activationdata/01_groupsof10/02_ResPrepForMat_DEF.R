# Structuur Matlab file
# For a single axial slice (Z = 2 mm) of data, you should have:
# Bmap_N_S: 'Difference between Novel and Standard betas averaged over 28 subjects'
# Tmap_N_S: 'T-statistics for the paired t-test comparing Novel and Standard betas'
# Pmap_N_S: 'Binary map indicating significance at P<0.001 (fdr corrected)'
# Underlay: 'Structural image ch2bet from MRIcron, warped to functional data'
rm(list=ls())
library(oro.nifti)
library(R.matlab)
f <- c("01_N10_smooth_EXCL/")#, "01_N10_smooth6_EXCL/", "02_N20_smooth_EXCL/", "02_N20_smooth6_EXCL/")
f <- c("01_N10_smooth6_EXCL/")
p <- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/01_smallN/01_results/"
s.p <- "./Onderzoek/doctoraat/reports_abstracts/Paper_finaal/03_manuscript/00_pix/"
ZZ=i=604
path1 <- paste(p,f, sep="")
path2 <- paste(path1, sprintf("%04d", ZZ), sep="")


s.path  <- "./Onderzoek/doctoraat/scripts/Paper04/01_data/00_singleSubject/02_results/"
s.path2 <- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_ALL/results/"
a.path <-  "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/01_analysis/MNI152_T1_2mm_brain.nii.gz"

#### LET OP FSL START OP VOXEL ID == 00 !!! C codering en NIET op 1!
#Low RES
i=43
B1 <- readNIfTI(paste(path2, "/cl.nii.gz", sep=""))
Z1 <- readNIfTI(paste(path2, "/zstat1.nii.gz", sep=""))
S1 <- readNIfTI(paste(path2, "/OUTVAR.nii.gz", sep=""))
A1 <- readNIfTI(a.path)[,,,1]

z <- Z1[,,i]
b <- B1[,,i]
a <- A1[,,i]
s <- S1[,,i]

filename <- paste(s.p, "CON_emoS6_", i-1,".mat", sep="")
writeMat(filename, A=a, B=b, S=s, Z=z)

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
