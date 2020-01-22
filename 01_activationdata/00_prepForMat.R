# Structuur Matlab file
# For a single axial slice (Z = 2 mm) of data, you should have:
# Bmap_N_S: 'Difference between Novel and Standard betas averaged over 28 subjects'
# Tmap_N_S: 'T-statistics for the paired t-test comparing Novel and Standard betas'
# Pmap_N_S: 'Binary map indicating significance at P<0.001 (fdr corrected)'
# Underlay: 'Structural image ch2bet from MRIcron, warped to functional data'
rm(list=ls())
library(oro.nifti)
library(R.matlab)
library(AnalyzeFMRI)

# zfile 		<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/00_singleSubject/01/tfMRI_EMOTION_LR/stats/zstat3.nii.gz"
# bfile 		<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/00_singleSubject/01/tfMRI_EMOTION_LR/stats/outSingSubj-con3.nii.gz"
# stab.file 	<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/00_singleSubject/01/out/OUT-con3.nii.gz"
# group.file	<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/01_analysis/results/OUT77_3.nii.gz"
# groupz.file	<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/01_analysis/results/zstat1.nii.gz"
# afile 		<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/00_singleSubject/01/out/MNI152_T1_2mm_brain.nii.gz"
 s.path 		<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/00_singleSubject/01/out/"
#### LET OP FSL START OP VOXEL ID == 00 !!! C codering en NIET op 1!
#Low RES
#system(paste("gunzip", zfile))
#system(paste("gunzip", bfile))
#system(paste("gunzip", stab.file))
#system(paste("gunzip", group.file))
#system(paste("gunzip", groupz.file))
#system(paste("gunzip", afile))
zfile 		<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/00_singleSubject/01/tfMRI_EMOTION_LR/stats/zstat3.nii"
bfile 		<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/00_singleSubject/01/tfMRI_EMOTION_LR/stats/outSingSubj-con3.nii"
stab.file 	<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/00_singleSubject/01/out/OUT-con3.nii"
group.file	<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/01_analysis/results/OUT77_3.nii"
groupz.file	<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/01_analysis/results/zstat1.nii"
afile 		<- "./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/00_singleSubject/01/out/MNI152_T1_2mm_brain.nii"
B1 <- readNIfTI(bfile)#, 		verbose=TRUE,warn=1,reorient=FALSE)
Z1 <- readNIfTI(zfile)#, 		verbose=TRUE,warn=1,reorient=FALSE)
S1 <- readNIfTI(stab.file)#,  verbose=TRUE,warn=1,reorient=FALSE)
A1 <- readNIfTI(afile)[,,,1]#, 		verbose=TRUE,warn=1,reorient=FALSE)[,,,1]
B2 <- readNIfTI(group.file)#, verbose=TRUE,warn=1,reorient=FALSE)
Z2 <- readNIfTI(groupz.file)#,verbose=TRUE,warn=1,reorient=FALSE)

for(i in 1:91){
	z  <- Z1[,,i]
	b  <- B1[,,i]
	a  <- A1[,,i]
	s  <- S1[,,i]
	b2 <- B2[,,i]
	z2 <- Z2[,,i]

	filename <- paste(s.path, "CON_emo_", i-1, ".mat", sep="")
	writeMat(filename, A=a, B=b, S=s, Z=z, B2=b2)
}
