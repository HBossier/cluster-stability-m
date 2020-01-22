#tmpdir='./Onderzoek/doctoraat/scripts/Paper04/01_data/1/MNINonLinear/Results/tfMRI_LANGUAGE_LR/tfMRI_LANGUAGE_LR_hp200_s4.feat/stats2++++/'
#tmpdir2='./Onderzoek/doctoraat/scripts/PaperReSelectionRate/01_data/1/MNINonLinear/Results/tfMRI_LANGUAGE_LR/tfMRI_LANGUAGE_LR_hp200_s4.feat/'
# set directories
tmpdir='./Onderzoek/Paper04/01_data/02_emotion'
savedir='./Onderzoek/Paper04/01_data/02_emotion/00_SingleSubject/'

# import neccessary libraries
import os as os
import sys 
import numpy as np
import nibabel as nib

# change directory and check whethere arguments are read in correctly
os.chdir(tmpdir)
os.listdir(tmpdir)
print(sys.argv)
# seed is passed through via the command line command and varies from 1:1000
seed=np.array(int(sys.argv[1])*1)
print(seed)
folder=str("run%04d/" % (seed))
print(folder)


### BOOTSTRAP
# bootstrap from unwhitened residuals
resids = os.path.join(tmpdir, 'res4d.nii.gz')
Res = nib.nifti1.load(resids)

# set random seed
np.random.seed(seed+1-1)
# choose random starting position with a block length of 16, 20 times (last block is only 12)
Tidx = np.random.choice(range(0,160,16), 11)
Tboot = [] 
for i in range(0,175,16): 
    j = i/16;
    Tboot[i:i+16] = range(Tidx[j],Tidx[j]+16,1)
    
# select only the 316 first elements 
#Tboot = Tboot[0:176]
# get the data effectively
data = Res.get_data()
data = data[:,:,:,Tboot]
affine = Res.get_affine()
#save the data
boot = nib.Nifti1Image(data, affine, Res.get_header())
#clean up memory
del data
print("data is shuffled")

### DESIGN MATRIX + BETAS
# read in the design matrix
X=np.loadtxt('design.mat', skiprows=5)
X=np.array(X)
# check whether both have same size
# print(X.shape)

#read in the beta images
BETA1 = nib.load(os.path.join(tmpdir, 'pe1.nii.gz'))
B1 = BETA1.get_data()
AFFINE= BETA1.get_affine()
b1 = B1.reshape((np.product(BETA1.shape)))
X1=X[:,0]
xb1=np.outer(b1,X1)
#XB1=xb1.reshape((91,109,91,316), order="F")
#xb1img = nib.Nifti1Image(XB1, AFFINE)
del B1


BETA2 = nib.load(os.path.join(tmpdir, 'pe2.nii.gz'))
B2 = BETA2.get_data()
b2 = B2.reshape((np.product(BETA2.shape)))
X2=X[:,1]
xb2=np.outer(b2,X2)
#XB2=xb2.reshape((91,109,91,316), order="F")
#xb2img = nib.Nifti1Image(XB2, AFFINE)
del B2


BETA3 = nib.load(os.path.join(tmpdir, 'pe3.nii.gz'))
B3 = BETA3.get_data()
b3 = B3.reshape((np.product(BETA3.shape)))
X3=X[:,1]
xb3=np.outer(b3,X3)
#XB3=xb3.reshape((91,109,91,316), order="F")
#xb3img = nib.Nifti1Image(XB3, AFFINE)
del B3

BETA4 = nib.load(os.path.join(tmpdir, 'pe4.nii.gz'))
B4 = BETA4.get_data()
b4 = B4.reshape((np.product(BETA4.shape)))
X4=X[:,1]
xb4=np.outer(b4,X4)
#XB4=xb4.reshape((91,109,91,316), order="F")
#xb4img = nib.Nifti1Image(XB4, AFFINE)
del B4

print("set design matrix")
# set the activation  (mean image is added in the bash script)
outtmp=xb1+xb2+xb3+xb4
print("design matrix is set.")
del xb1, xb2, xb3, xb4
print("load bootstrap data")
outtmp2=boot.get_data()
# for computational reasons data is reshaped
out= outtmp+outtmp2.reshape((902629,176))
del outtmp, outtmp2
print("all set for saving.")

# and than reshaped again.
out=out.reshape((91,109,91,176))
new_img = nib.Nifti1Image(out,AFFINE,Res.get_header())
print("start saving file ..")
nib.save(new_img, savedir+folder+'new_img2.nii.gz')
print("\n .")

