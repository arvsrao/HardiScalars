#!/usr/bin/env python

import sys, os, nifti, math
from numpy import *

if len(sys.argv) < 2: 
   print "Usage:",sys.argv[0],"subject #","\n"
   sys.exit(1)
else:   
   id = sys.argv[1]
   DIR = sys.argv[2]
   TPLfile = sys.argv[3] # template first peak file. 
   
first_peak_subj=nifti.NiftiImage( DIR + '/' + id + '0.nii.gz')
second_peak_subj=nifti.NiftiImage( DIR + '/' + id + '1.nii.gz')

first_peak_tple=nifti.NiftiImage( TPLfile )
#second_peak_tple=nifti.NiftiImage('/sbia/home/raoarvin/comp_space/autism_maps/R0090/Norm_FOD/R0090_Norm_FOD_Warp_first_peak.nii.gz')

wm=nifti.NiftiImage('/sbia/home/raoarvin/comp_space/autism_maps/R0090/Norm_FOD/R0090_wmMask.nii.gz')
narray1=zeros((3,1, 80,128,128),float32)
narray2=zeros((3,1, 80,128,128),float32)

a=first_peak_subj.data[:,0,:,:,:]
b=second_peak_subj.data[:,0,:,:,:]
c=first_peak_tple.data[:,0,:,:,:]

norm_a_map=sqrt( a[0,:,:,:]*a[0,:,:,:] + a[1,:,:,:]*a[1,:,:,:] + a[2,:,:,:]*a[2,:,:,:] )
norm_b_map=sqrt( b[0,:,:,:]*b[0,:,:,:] + b[1,:,:,:]*b[1,:,:,:] + b[2,:,:,:]*b[2,:,:,:] )
norm_c_map=sqrt( c[0,:,:,:]*c[0,:,:,:] + c[1,:,:,:]*c[1,:,:,:] + c[2,:,:,:]*c[2,:,:,:] )

for i in range(0,80):
    for j in range(0,128):
        for k in range(0,128):
            if wm.data[i,j,k]!=0:
                 dot_a = dot( a[:,i,j,k]/norm_a_map[i,j,k], c[:,i,j,k]/norm_c_map[i,j,k] )
                 dot_b = dot( b[:,i,j,k]/norm_b_map[i,j,k], c[:,i,j,k]/norm_c_map[i,j,k] )
                 if dot_b > dot_a and norm_b_map[i,j,k] > .1: #when do you have a second peak (.10 threshold ) and is it aligned with the first peak of the template.
                     narray2[:,0,i,j,k]= first_peak_subj.data[:,0,i,j,k]
                     narray1[:,0,i,j,k]= second_peak_subj.data[:,0,i,j,k]
                 else:
                     narray1[:,0,i,j,k]= first_peak_subj.data[:,0,i,j,k]
                     narray2[:,0,i,j,k]= second_peak_subj.data[:,0,i,j,k]
            else:
                 narray1[:,0,i,j,k]=array([0,0,0], dtype=float32)
                 narray2[:,0,i,j,k]=array([0,0,0], dtype=float32)

h=first_peak_tple.asDict()

nim1=nifti.NiftiImage(narray1*wm.data, header=h)
nim1.save( DIR + '/' + id + '_relabeled_1st_peak.nii.gz')

nim2=nifti.NiftiImage(narray2*wm.data, header=h)
nim2.save( DIR + '/' + id + '_relabeled_2nd_peak.nii.gz')
