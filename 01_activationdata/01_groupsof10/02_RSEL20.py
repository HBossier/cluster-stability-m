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
seed=np.array(tmp1)
maxval=int(sys.argv[2]) * 1

# smooth sigma in mm 1.69851380042
np.random.seed(seed+1-1)
tmp=np.random.choice(range(0,maxval,1), 20, replace=0)
print ' '.join(map(str, tmp))
