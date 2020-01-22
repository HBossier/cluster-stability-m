scp ***./Onderzoek/Paper04/02_results/00_ResRates/01_smallN/01_Noise/OUT* .

fslmerge -t copeREG `
for i in 11 12 13 14 15
do
echo ./$i/cope4.nii
done
`
# merge the variance files
fslmerge -t varcopeREG `
for i in 11 12 13 14 15
do
echo $i/varcope4.nii
done
`

# merge the masks
fslmerge -t MASK `
for i in 11 12 13 14 15
do
echo ./$i/mask.nii
done
`
$FSLDIR/bin/flameo --cope=./copeREG --varcope=./varcopeREG --mask=MASK --dm=design.mat --tc=design.con --cs=design.grp --runmode=flame1 --ld=./01_OrRes
$FSLDIR/bin/cluster -i ./01_OrRes/zstat1.nii.gz -t 3.090232 --oindex=./01_OrRes/cl-con > ./01_OrRes/con


foldid='./Onderzoek/doctoraat/scripts/Paper04/01_data/01_smallN/01_origData/01_OrRes'
$FSLDIR/bin/fslmaths $foldid/zstat1.nii.gz  -mas mask -cpval -g smooothness pmapFWE
$FSLDIR/bin/fdr -i pmap -m mask -a FDRcorrected


fslmerge -t copeREG `
for i in 11 12 13 14 15
do
echo /Volumes/SANNE_HD/hcp_data/$i/MNINonLinear/Results/tfMRI_LANGUAGE_LR/tfMRI_LANGUAGE_LR_hp200_s4.feat/stats/cope3.nii.gz
done
`

# merge the variance files
fslmerge -t varcopeREG `
for i in 11 12 13 14 15
do
echo /Volumes/SANNE_HD/hcp_data/$i/MNINonLinear/Results/tfMRI_LANGUAGE_LR/tfMRI_LANGUAGE_LR_hp200_s4.feat/stats/varcope3.nii.gz
done
`
fslmerge -t MASK `
for i in 11 12 13 14 15
do
echo /Volumes/SANNE_HD/hcp_data/$i/MNINonLinear/Results/tfMRI_LANGUAGE_LR/tfMRI_LANGUAGE_LR_hp200_s4.feat/mask.nii.gz
done
`


$FSLDIR/bin/flameo --cope=./copeREG --varcope=./varcopeREG --mask=MASK --dm=../../design.mat --tc=../../design.con --cs=../../design.grp --runmode=flame1 --ld=./results
$FSLDIR/bin/cluster -i ./results/zstat1.nii.gz -t 3.090232 --oindex=./results/cl-con > ./results/con
