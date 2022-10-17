#!/usr/bin/env python

import numpy as np
import math
import argparse

parser = argparse.ArgumentParser(description="derive spectral thermal conductivity of coherent component")

parser.add_argument("fin",type=str,help="kc_elem file from alamode at a temperature")
parser.add_argument("fout",type=str,default="spec.dat",help="spectral thermal conductivity data")
parser.add_argument('-nb','--nbranches',type=int,help="the number of phonon branches")
parser.add_argument('-nq','--nqpoints',type=int,help='the number of irreducible q points')
args = parser.parse_args() 

nb=args.nbranches
nq=args.nqpoints
spec=np.zeros([nq,nb,4])
with open(args.fin) as f:
    for line in f:
        element=line.strip().split()
        ib=int(element[5])
        iq=int(element[3])
        ix=int(element[1])
        spec[ib,iq,0] = float(element[6])
        spec[ib,iq,ix] += float(element[8])

spec[:,:,0]=spec[:,:,0]/33.356
spec.reshape([nb*nq,4])
index=np.argsort(spec[:,0])
spec_sorted=spec[index,:]
freq_max=spec_sorted[-1,0]
df=0.1
nf=math.ceil(freq_max/df)
spec_kappa=np.zeros([nf,4])
ifreq=1
for i in range(len(spec_sorted)):
    if spec_sorted[i,0] <=  ifreq*0.1:
        spec_kappa[ifreq-1,:] += spec_sorted[i,:]
    else:
        ifreq += 1
        spec_kappa[ifreq-1,:] += spec_sorted[i,:]
        
accu_kappa=np.zeros([nf,4])
accu_kappa[0,:]=spec_sorted[0,:]
for i in range(len(spec_sorted-1)):
    accu_kappa[i+1,0]=spec_sorted[i+1,0]
    accu_kappa[i+1,1:4]=accu_kappa[i,1:4]+spec_sorted[i+1,1:4]

with open(args.fout,w) as f:
    for i in range(nq*nb):
        f.write("%15.5e %15.5e %15.5e %15.5e " % tuple(spec_sorted[i,:]))
        f.write("%15.5e %15.5e %15.5e\n" % tuple(accu_kappa[i,1:4]))






