#corrected p values FWE
path1="./Onderzoek/doctoraat/scripts/Paper04/01_data/02_emotion/01_smallN/01_results/00"
for i in {1..50}
do
	echo $i 
	simid=$(printf '%02d' $i)
	cd $path1$simid
	awk '$1~/VOLUME/{v=$2};$1~/RESELS/{r=$2};END{printf("%g",1.0*v/r)}' smoothness00$simid > tmp
	RESELcount=$(cat tmp)
	thres=$(ptoz 0.2 -g $RESELcount)
	fslmaths zstat1.nii.gz -thr $thres thresholded.nii.gz
done

awk '$1~/VOLUME/{v=$2};$1~/RESELS/{r=$2};END{printf("%g",1.0*v/r)}' smoothness > tmp
RESELcount=$(cat tmp)
thres=$(ptoz 0.05 -g $RESELcount)
fslmaths zstat1.nii.gz -thr $thres thresholded.nii.gz


