for (( i=1; i<=6; i++ ))
	do	
	flirt -in ./Onderzoek/doctoraat/scripts/Paper04/01_data/1/out/out${i}.nii.gz -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -out ./Onderzoek/doctoraat/scripts/Paper04/01_data/02_results/outBIN${i}HIGH -omat ./Onderzoek/doctoraat/scripts/Paper04/01_data/02_results/outbin${i}.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear
	done
