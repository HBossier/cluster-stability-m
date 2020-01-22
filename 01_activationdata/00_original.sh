film_gls --rn=stats2 --sa --ms=5 --in=filtered_func_data --pd=design.mat 1000.0 -v


film_gls --rn=stats --in=tfMRI_EMOTION_LR --thr=1000.0 --con=tfMRI_EMOTION_LR_hp200_s4_level1.con --pd=tfMRI_EMOTION_LR_hp200_s4_level1.mat -v 


$FSLDIR/bin/cluster -i zstat3.nii.gz -t 3.090232 --oindex=CL_INDEX_c3 > CL_c3
Rscript ../../../../../../02_scripts/01_language/cluster.R CL_c3
min=$(awk 'FNR == 1 {print $1}' tst)
fslmaths CL_INDEX_c3 -thr $min -bin outSingSubj-con3

fslmaths CL_INDEX_c3 -thr $min outSingSubj-con3-values


fslmaths tfMRI_EMOTION_LR -Tmean mean_func_filtered