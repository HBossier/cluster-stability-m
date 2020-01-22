# tmpdir='./Onderzoek/Paper04/01_data'
# savedir='./Onderzoek/Paper04/01_data/'
# # import neccessary libraries
import os as os
import sys
import numpy as np
import math

# # seed is passed through via the command line command and varies from 1:1000
#print(sys.argv)
tmp1=int(sys.argv[1]) * 1
seed=np.array(tmp1)//2 
maxval=int(sys.argv[2]) * 1

np.random.seed(seed+1-1)
tmp=np.random.choice(maxval, 40, replace=0)
if tmp1%2==0 :
    print ' '.join(map(str, tmp[0:20]))
else:
    print ' '.join(map(str, tmp[20:40]))