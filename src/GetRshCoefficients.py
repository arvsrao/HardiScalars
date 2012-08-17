#!/usr/bin/env python

#
#  Extracts the RSH coefficents for a particular voxel from a given RSH image.
# 
#  User must input the voxel address of interest. 
#

                
import sys,os
import nifti
import numpy

if len(sys.argv) < 2: 
   print "Usage:",sys.argv[0],"pixel coordinates: x, y, z  /   *.hdr (or nii.gz)\n"
   sys.exit(1)

X = sys.argv[1]
Y = sys.argv[2]
Z = sys.argv[3]
      
for i in sys.argv[4:]:
  if not os.path.exists(i):
      print i,"does not exist!"
      sys.exit(1)
  im=nifti.NiftiImage(i)

#numpy.set_printoptions(precision=10)
#numpy.set_printoptions(suppress=True)
print 'Coefficients are: \n'
try:
   print [j for j in im.data[:,0,Z,Y,X]] 
except IndexError, e:
   print "Data out of range for", i, "\n\nThe error is:",e
   

sys.exit(0)
