# cluster -based with 75 cut-off

## no noise
foldid2="./Onderzoek/Paper04/02_results/02_smallN/00_noNoise"
for (( i=1; i<=1000; i++ ))	
do
foldid=$(printf './Onderzoek/Paper04/01_data/01_smallN/00_data/run%04d' $i)
cluster -i $foldid/00_noNoise/zstat1.nii.gz -t 3.090232 --oindex=$foldid2/cl$i-con$j > $foldid2/f$i-con$j
done

#noise
foldid2="./Onderzoek/Paper04/02_results/02_smallN/01_Noise"
for (( i=1; i<=1000; i++ ))	
do
foldid=$(printf './Onderzoek/Paper04/01_data/01_smallN/00_data/run%04d' $i)
cluster -i $foldid/01_Noise/zstat1.nii.gz -t 3.090232 --oindex=$foldid2/cl$i-con$j > $foldid2/f$i-con$j
done


# give id for folder
## no noise
foldid3="./Onderzoek/Paper04/02_results/00_ResRates/01_smallN/00_noNoise"
for (( i=1; i<=1000; i++ ))	
do
foldid2="./Onderzoek/Paper04/02_results/02_smallN/00_noNoise"
Rscript cluster.R $foldid2/f$i-con
min=$(awk -v j=$j 'FNR == j {print $1}' tst)
fslmaths $foldid2/cl$i-con.nii.gz -thr $min -bin $foldid3/out$i-con
done

## noise
foldid3="./Onderzoek/Paper04/02_results/00_ResRates/01_smallN/01_Noise"
for (( i=1; i<=1000; i++ ))	
do
foldid2="./Onderzoek/Paper04/02_results/02_smallN/01_Noise"
Rscript cluster.R $foldid2/f$i-con
min=$(awk -v j=$j 'FNR == j {print $1}' tst)
fslmaths $foldid2/cl$i-con.nii.gz -thr $min -bin $foldid3/out$i-con
done


# make empty volumes with a little trick 
i=1
foldid3="./Onderzoek/Paper04/02_results/00_ResRates/01_smallN/00_noNoise"
fslmaths $foldid3/out$i-con -thr 10 -bin $foldid3/OUT-con
foldid4="./Onderzoek/Paper04/02_results/00_ResRates/01_smallN/01_Noise"
fslmaths $foldid4/out$i-con -thr 10 -bin $foldid4/OUT-con


for (( i=1; i<=1000; i++ ))	
do
echo $i
fslmaths $foldid3/OUT-con -add $foldid3/out$i-con $foldid3/OUT-con
fslmaths $foldid4/OUT-con -add $foldid4/out$i-con $foldid4/OUT-con
done


$FSLDIR/bin/fslmaths zstat1.nii.gz -ztop p1WLS
# make the files to save the number of selected voxels and the number of CORRECTLY selected voxels.
touch FDR_WLS_all.txt

## no noise
odir='./Onderzoek/Paper04/02_scripts'
foldid2="./Onderzoek/Paper04/02_results/02_smallN/00_noNoise"
for (( i=1; i<=1000; i++ ))	
do
	echo $i
	foldid=$(printf './Onderzoek/Paper04/01_data/01_smallN/00_data/run%04d' $i)
	$FSLDIR/bin/fslmaths $foldid/00_noNoise/zstat1.nii.gz -ztop $foldid2/pFDR$i
	$FSLDIR/bin/fdr -i $foldid2/pFDR$i -m $odir/maks_orig -q 0.1 > FDRtresh
	FDRt=$(echo $(cat FDRtresh) | awk '{ v=$4 };END{printf("%g", 1.0-1.0*v)}')
	# obtain a comparable map with (higher = significant inversion of the p-map via)
	$FSLDIR/bin/fslmaths $foldid2/pFDR$i -mul -1 -add 1 -thr $FDRt -mas  $odir/maks_orig -bin $foldid2/FDR$i.nii.gz
done

## no noise
odir='./Onderzoek/Paper04/02_scripts'
foldid2="./Onderzoek/Paper04/02_results/02_smallN/01_Noise"
for (( i=1; i<=1000; i++ ))	
do
	echo $i
	foldid=$(printf './Onderzoek/Paper04/01_data/01_smallN/00_data/run%04d' $i)
	$FSLDIR/bin/fslmaths $foldid/01_Noise/zstat1.nii.gz -ztop $foldid2/pFDR$i
	$FSLDIR/bin/fdr -i $foldid2/pFDR$i -m $odir/maks_orig -q 0.1 > FDRtresh
	FDRt=$(echo $(cat FDRtresh) | awk '{ v=$4 };END{printf("%g", 1.0-1.0*v)}')
	# obtain a comparable map with (higher = significant inversion of the p-map via)
	$FSLDIR/bin/fslmaths $foldid2/pFDR$i -mul -1 -add 1 -thr $FDRt -mas $odir/maks_orig -bin $foldid2/FDR$i.nii.gz
done



for i in $(cat alpha)
do
	# computes the threshold
	$FSLDIR/bin/fdr -i p1WLS -m ./mask -q $i > FDRtresh
	#stores the threshold 1-P !!!
	FDRt=$(echo $(cat FDRtresh) | awk '{ v=$4 };END{printf("%g", 1.0-1.0*v)}')
	# obtain a comparable map with (higher = significant inversion of the p-map via)
	base=${i##*.}
	fname=$(echo FDR_1minp$base)
	$FSLDIR/bin/fslmaths p1WLS -mul -1 -add 1 -thr $FDRt -mas ./mask ${fname}.nii.gz
	# determine the number and save it (ALL detected voxels)
	tmp=$(fslstats ${fname} -V)
	echo $tmp | cat - FDR_WLS_all.txt > temp && mv temp FDR_WLS_all.txt
	#clean -up
	rm ${fname}.nii.gz
done
