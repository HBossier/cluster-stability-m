# first level more smoothing.
sed -i.orig -e 's/_s4/_s4_6/g' -e 's/fmri(smooth)\ 4/fmri(smooth)\ 6/g' tfMRI_EMOTION_LR_hp200_s4_level1.fsf
mv tfMRI_EMOTION_LR_hp200_s4_level1.fsf  tfMRI_EMOTION_LR_hp200_s4_level1_6.fsf 
mv tfMRI_EMOTION_LR_hp200_s4_level1.fsf.orig  tfMRI_EMOTION_LR_hp200_s4_level1.fsf 

# first level analysis.
path1='./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/00_rawData/'
cd $path1
files=(`cat "files.txt"`) 
i=1
for i in {2..76}
do
	echo $i 
	j=$(expr $i + 1)
	cd $path1$j/${files[$i]}/MNINonLinear/Results/tfMRI_EMOTION_LR
	sed -i.orig -e 's/_s4/_s4_6/g' -e 's/fmri(smooth)\ 4/fmri(smooth)\ 6/g' tfMRI_EMOTION_LR_hp200_s4_level1.fsf
	mv tfMRI_EMOTION_LR_hp200_s4_level1.fsf  tfMRI_EMOTION_LR_hp200_s4_level1_6.fsf 
	mv tfMRI_EMOTION_LR_hp200_s4_level1.fsf.orig  tfMRI_EMOTION_LR_hp200_s4_level1.fsf 
	feat_model tfMRI_EMOTION_LR_hp200_s4_level1_6
	feat tfMRI_EMOTION_LR_hp200_s4_level1_6.fsf
done

for i in {12..76}
do
	echo $i 
	j=$(expr $i + 1)
	cd $path1$j/${files[$i]}/MNINonLinear/Results/tfMRI_EMOTION_LR
	path2=$path1$j/${files[$i]}/MNINonLinear/Results/tfMRI_EMOTION_LR
	_r1="${path2//\//\\/}"
	sed -i -e "s/\.\//{$_r1}\//g" tfMRI_EMOTION_LR_hp200_s4_level1.fsf
	sed -i -e "s/{//g" tfMRI_EMOTION_LR_hp200_s4_level1.fsf
	sed -i -e "s/}//g" tfMRI_EMOTION_LR_hp200_s4_level1.fsf
	feat_model tfMRI_EMOTION_LR_hp200_s4_level1	
	feat tfMRI_EMOTION_LR_hp200_s4_level1
done

 
## Escape path for sed using bash find and replace 
_r1="${_r1//\//\\/}"

# first level with UNSMOOTHED DATA
path1='./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/00_rawData/'
cd $path1
files=(`cat "files.txt"`) 
for i in {0..76}
do
	j=$(expr $i + 1)
	cd $path1$j/${files[$i]}/MNINonLinear/Results/tfMRI_EMOTION_LR
	sed -i-e 's/\.\./\./g' tfMRI_EMOTION_LR_hp200_s4_level1.fsf
	feat_model tfMRI_EMOTION_LR_hp200_s4_level1	
	film_gls --rn=stats --in=tfMRI_EMOTION_LR --thr=1000.0 --con=tfMRI_EMOTION_LR_hp200_s4_level1.con --pd=tfMRI_EMOTION_LR_hp200_s4_level1.mat -v
done


# SECOND level analysis WITHOUT SMOOTHING
path1='./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/'
cd $path1
files=(`cat "files.txt"`) 
fslmerge -t copeREG `
for i in {0..76}
do
	j=$(expr $i + 1)
	echo $j/${files[$i]}/MNINonLinear/Results/tfMRI_EMOTION_LR/stats/cope3
done
`

# merge the variance files
fslmerge -t varcopeREG `
for i in {0..76}
do
	j=$(expr $i + 1)
	echo $j/${files[$i]}/MNINonLinear/Results/tfMRI_EMOTION_LR/stats/varcope3
done
`
fslmerge -t MASK `
for i in {0..76}
do
	j=$(expr $i + 1)
	echo $j/${files[$i]}/MNINonLinear/Results/tfMRI_EMOTION_LR/*mask*
done
`


$FSLDIR/bin/flameo --cope=./copeREG --varcope=./varcopeREG --mask=MASK --dm=design.mat --tc=design.con --cs=design.grp --runmode=flame1 --ld=./results
$FSLDIR/bin/cluster -i ./results/zstat1.nii.gz -t 3.090232 --oindex=CL_3 > f_CL_3

Rscript ../../../../../02_scripts/02_emotion/cluster.R f_CL_3 $(pwd)
min=$(awk -v j=1 'FNR == j {print $1}' tst)
fslmaths CL_3.nii.gz -thr $min -bin OUT77_3


# SECOND LEVEL analysis with smoothing.
path1='./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/00_rawData/'
cd $path1
files=(`cat "files.txt"`) 
cd ./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/01_analysis/01_SMOOTH
fslmerge -t copeREG `
for i in {0..76}
do
	j=$(expr $i + 1)
	echo $path1$j/${files[$i]}/MNINonLinear/Results/tfMRI_EMOTION_LR/tfMRI_EMOTION_LR_hp200_s4.feat/stats/cope3
done
`

# merge the variance files
fslmerge -t varcopeREG `
for i in {0..76}
do
j=$(expr $i + 1)
echo $path1$j/${files[$i]}/MNINonLinear/Results/tfMRI_EMOTION_LR/tfMRI_EMOTION_LR_hp200_s4.feat/stats/varcope3
done
`

fslmerge -t MASK `
for i in {0..76}
do
j=$(expr $i + 1)
echo $path1$j/${files[$i]}/MNINonLinear/Results/tfMRI_EMOTION_LR/tfMRI_EMOTION_LR_hp200_s4.feat/*mask*
done
`
$FSLDIR/bin/flameo --cope=./copeREG --varcope=./varcopeREG --mask=MASK --dm=design.mat --tc=design.con --cs=design.grp --runmode=flame1 --ld=./results
$FSLDIR/bin/cluster -i ./results/zstat1.nii.gz -t 3.090232 --oindex=CL_3 > f_CL_3

Rscript ../../../../../02_scripts/02_emotion/cluster.R f_CL_3 $(pwd)
min=$(awk -v j=1 'FNR == j {print $1}' tst)
fslmaths CL_3.nii.gz -thr $min -bin OUT75_3


# compute the intrinsic smoothness of the HCP data based on the original residuals (SMOOTHED VERSION, unsmoothed => remove tfMRI... from folder path.). 
path1='./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/00_rawData/'
cd $path1
files=(`cat "files.txt"`) 
i=1
for i in {0..76}
do
echo $i 
j=$(expr $i + 1)
cd $path1$j/${files[$i]}/MNINonLinear/Results/tfMRI_EMOTION_LR/tfMRI_EMOTION_LR_hp200_s4.feat/stats
smoothest -d 172 -r res4d.nii.gz -m ../*mask* -V > smoothness_SM
done

# compute the intrinsic smoothness of the HCP data based on the original residuals (SMOOTHED VERSION, unsmoothed => remove tfMRI... from folder path.). 
path1='./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/02_all/00_rawData/'
cd $path1
files=(`cat "files.txt"`) 
i=1
for i in {0..76}
do
echo $i 
j=$(expr $i + 1)
cd $path1$j/${files[$i]}/MNINonLinear/Results/tfMRI_EMOTION_LR/tfMRI_EMOTION_LR_hp200_s4.feat/stats
smoothest -d 172 -r res4d.nii.gz -m ../*mask* -V > smoothness_SM
done