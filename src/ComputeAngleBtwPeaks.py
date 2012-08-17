#!/usr/bin/env python

import sys, os, nifti, math
from numpy import *

if len(sys.argv) < 2: 
   print "Usage:",sys.argv[0],"subject #","\n"
   sys.exit(1)
else:   
   id = sys.argv[1]

first_peak_subj=nifti.NiftiImage('/sbia/home/raoarvin/comp_space/GeoVariance_Maps/'+id+'/Norm_FOD/'+id+'_
.nii.gz')
second_peak_subj=nifti.NiftiImage('/sbia/home/raoarvin/comp_space/GeoVariance_Maps/'+id+'/Norm_FOD/'+id+'_second_peak.nii.gz')

first_peak_tple=nifti.NiftiImage('/sbia/home/raoarvin/comp_space/GeoVariance_Maps/R0090/Norm_FOD/R0090_first_peak.nii.gz')
second_peak_tple=nifti.NiftiImage('/sbia/home/raoarvin/comp_space/GeoVariance_Maps/R0090/Norm_FOD/R0090_second_peak.nii.gz')

wm=nifti.NiftiImage('/sbia/home/raoarvin/comp_space/GeoVariance_Maps/R0090/Norm_FOD/R0090_wmMask.nii.gz')
narray=zeros((80,128,128),float32)

a=first_peak_subj.data[:,0,:,:,:]
b=second_peak_subj.data[:,0,:,:,:]
c=first_peak_tple.data[:,0,:,:,:]

norm_a_map=sqrt( a[0,:,:,:]*a[0,:,:,:] + a[1,:,:,:]*a[1,:,:,:] + a[2,:,:,:]*a[2,:,:,:])
norm_b_map=sqrt( b[0,:,:,:]*b[0,:,:,:] + b[1,:,:,:]*b[1,:,:,:] + b[2,:,:,:]*b[2,:,:,:])
norm_c_map=sqrt( c[0,:,:,:]*c[0,:,:,:] + c[1,:,:,:]*c[1,:,:,:] + c[2,:,:,:]*c[2,:,:,:])

for i in range(0,80):
    for j in range(0,128):
        for k in range(0,128):
            if wm.data[i,j,k]!=0:
                 dot_a = dot( a[:,i,j,k]/norm_a_map[i,j,k], c[:,i,j,k]/norm_c_map[i,j,k] )
                 dot_b = dot( b[:,i,j,k]/norm_b_map[i,j,k], c[:,i,j,k]/norm_c_map[i,j,k] )
                 if dot_a > dot_b:
                     narray[i,j,k]=math.degrees(math.acos(dot_a))
                 else:
                     narray[i,j,k]=math.degrees(math.acos(dot_b))
            else:
                 narray[i,j,k]=0      

h=first_peak_tple.asDict()
h['dim'][5]=1
h['intent_code']=0

#nim1=nifti.NiftiImage(narray, header=h)
#nim1.save('/sbia/home/raoarvin/comp_space/GeoVariance_Maps/'+id+'/Norm_FOD/'+id+'_orientation_diff_masked.nii.gz')

nim2=nifti.NiftiImage(narray*wm.data, header=h)
nim2.save('/sbia/home/raoarvin/comp_space/GeoVariance_Maps/'+id+'/Norm_FOD/'+id+'_first_peak_orientation_diff_masked.nii.gz')
