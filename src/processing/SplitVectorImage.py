#!/usr/bin/env python

import sys, os, nifti

if len(sys.argv) < 3: 
   print "Usage:",sys.argv[0],"vector image *.hdr or .nii \  prefix filenames \n"
   sys.exit(1)
   
prefix = sys.argv[2]
print prefix

for i in sys.argv[1:]:
  if os.path.exists(i):
      print i,"does exist!"
      im=nifti.NiftiImage(i)
      break
  else:
      print i, "doesn't exist"
      sys.exit(1)

h=im.asDict()
h['dim'][5]=1

nim1=nifti.NiftiImage(im.data[0,0,:,:,:], header=h)
nim2=nifti.NiftiImage(im.data[1,0,:,:,:], header=h)
nim3=nifti.NiftiImage(im.data[2,0,:,:,:], header=h)
nim4=nifti.NiftiImage(im.data[3,0,:,:,:], header=h)

nim1.save( prefix + '_first_peak_GC.nii')          
nim2.save( prefix + '_second_peak_GC.nii')
nim3.save( prefix + '_first_peak_mag.nii')
nim4.save( prefix + '_second_peak_mag.nii')
