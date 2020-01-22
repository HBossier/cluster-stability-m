# # import neccessary libraries
import os as os
import sys
import numpy as np
import math

# # seed is passed through via the command line command and varies from 1:1000
#print(sys.argv)
selection =  np.loadtxt(sys.argv[1], dtype=int)

tmp1=int(sys.argv[2]) * 1
tmp2=int(sys.argv[3]) * 1 
tmp3=int(sys.argv[4]) * 1 
seed=np.array(tmp1 + 500 * (tmp2-1))

# smooth sigma in mm 1.69851380042
np.random.seed(seed+1-1)
tmp=np.random.choice(range(0,tmp3,1), tmp3)
print ' '.join(map(str, selection[tmp]))
