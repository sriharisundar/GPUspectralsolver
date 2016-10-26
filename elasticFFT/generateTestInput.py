import numpy as np
import sys

filename=sys.argv[1];

dim=np.zeros(3,dtype=int)
for i in range(3):
	dim[i]=sys.argv[(2+i)]

eul=np.zeros(3)
for i in range(3):
	eul[i]=sys.argv[(5+i)]

output=open(filename,'w')

prodDim=np.prod(dim)

for k in range(dim[2]):
	for j in range(dim[1]):
		for i in range(dim[0]):
			output.write("%3.3f %3.3f %f %d %d %d 1 1 \n" % (eul[0],eul[1],eul[2],i+1,j+1,k+1))