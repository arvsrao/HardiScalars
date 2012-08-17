#!/usr/bin/env python
#------------------------------------------------------------------------------#
#ConvertToMrtrixBasis.py
#Created by Harini Evani, and later updated by Arvind Rao on 5/08/2012
#Usage: ConvertToMrtrixBasis.py --in=<SBIA FOD File> --prefix=<outfile> [-h] [-v] 
#	Options: --in=       --- input SBIA RSH file ,expects an argument
#	         --prefix=   --- prefix for output file name
#		 -h         --- displays a help message
#	         -v         --- verbose option
#		 -V 	    --- Version information			
#------------------------------------------------------------------------------#
#This script converts a SBIA FOD file into MRTRIX basis format 
#------------------------------------------------------------------------------#
#Dependencies: Python, Numpy, SbiaUtils
#Version 2.0
#Last Update 05/08/2012
#------------------------------------------------------------------------------#

#Dependent Modules :
import sys, getopt, string, signal, struct, tempfile, os,os.path
from SbiaUtils import *
from math import sqrt
import numpy
  
#-------------------------------

version_number = "2.0" # - MODIFY VERSION HERE!
last_modified = "05-08-2012" # - MODIFY DATE HERE! 

#-------------------------------
#function definition for version
def version():
	print '''
		This is version ''' + version_number + ''' of the script ConvertToMrtrixBasis.py.
		This was last modified on ''' + last_modified
	sys.exit(0)

#-------------------------------
#function definition for help message
def help_message():
    print '''
 Usage: ConvertToMrtrixBasis.py --in=<SBIA FOD File> --prefix=<outfile> [-h] [-v] 
	Options: --in=        --- input SBIA RSH rep file      ( expected )
	         --prefix=    --- prefix for output file name  ( expected )
	         --outputDir  --- output directory             ( defaults to current directory )
		  -h          --- displays a help message
	          -v          --- verbose                      
		  -V          --- Version information '''	     
    sys.exit(0)

#-------------------------------
#function that creates a temp directory
def make_temp(temp_suffix, temp_prefix, temp_path, temp_dir=[]):
	try:
   		tdir_path = tempfile.mkdtemp(suffix=temp_suffix, prefix=temp_prefix, dir=temp_path)
	except OSError:
   		 print 'Error:  Couldn\'t create temp file.'
   		 sys.exit(0)
	return tdir_path

#-------------------------------
#function that removes a temp(?) directory
def remove_temp(temp_dir):
	try:
		shutil.rmtree(temp_dir)
	except OSError:
   		 print 'Error:  Couldn\'t remove temp directory.'
   		 sys.exit(0)

#-------------------------------
#function that parses input arguments
def parse(infile=[], outfile=[], verbose=[]):
         

	try:
               options, xarguments = getopt.getopt(sys.argv[1:], 'hvV', ['in=', 'prefix=', 'outputDir='])
	except getopt.GetoptError, err:
 	       print str(err)
	       help_message()
	
	if len(options) < 1:
	   print '\n\t>> You are missing the input file name and prefix for the output!!\n'
	   sys.exit(0)  	

        #parse short options
	for a in options[:]:
    		if a[0] == '-h':
       		 	help_message()

        verbose=False
	for a in options[:]:
    		if a[0] == '-v':
        		print 'output to the screen will be verbose'
			verbose = True
        		options.remove(a)
        		
	for a in options[:]:
    		if a[0] == '-V':
        		version()

       #parse long options
	for a in options[:]:
    		if a[0] == '--in' and a[1] != '':
			infile = a[1] 
        		options.remove(a)
        		break
    		elif a[0] == '--in' and a[1] == '':
        		print 'Please provide an input argument. Ex. ID**.nii.gz '
        		sys.exit(0)

        outDir = os.path.curdir + os.path.sep  #default output directory is the current directory.
        for a in options[:]:
              if a[0] == '--outputDir' and a[1] != '':
                 if os.path.exists(a[1]):
                    outDir = a[1]
                 elif a[0] == '--outputDir' and a[1] == '':
                      print 'If you intended to specify an output directory, please do so next time.'
                      sys.exit(0)
                 else:
                    print 'Your output directory does not exist. Using the current directory instead.'
                    outDir = os.path.curdir + os.path.sep
         
	for a in options[:]:
    		if a[0] == '--prefix' and a[1] != '':
			outfile = a[1] +'_mrtrix.nii.gz'
        		options.remove(a)
        		break
    		elif a[0] == '--prefix' and a[1] == '':
        		print 'Please provide a prefix for the output filename'
        		sys.exit(0)
        		
        		
        		
	#Get the absolute path name and base name for input and output
	#expand the files into absolute paths and check
	infile = os.path.realpath(infile)                
        outDir = os.path.realpath(outDir)
        outfile = outDir + '/' + outfile
	 
	if not os.path.isfile(infile):
	   print "*** File does not exist: " + infile
	   sys.exit(2)

        #check Nifti
	if not isNiftiFile(infile):
	   print "*** Subject image must be Nifti file!"
	   sys.exit(2)
	    
	# check output dir
	if not os.path.exists(outDir):
	   print outDir, "does not exist, creating it..."
	   makeDir(outDir)
	if not isWritable(outDir):
	   cryandexit("Output directory is not writable!", outDir) 		
       		

	return infile,outfile,verbose
#-------------------------------
#function that hooks the SIGINT interrupt(for Ctrl+C) to the signal handler(exit script)
def signal_handler(signal, frame):
        print 'You pressed Ctrl+C!'
        sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)
#-------------------------------

def main():
    
    mrtrix = True
    dataFile,outFile,verbose = parse() # parses argument and returns relavent values
    
    #print some info
    if verbose:
        print "dataFile: ", dataFile

    #check image header stuff
    data_im = NiftiImage(dataFile)
    data_intent = data_im.raw_nimg.intent_code

    try:
        if mrtrix:
           SbiaFodToMRtrixRun(data_im,outFile,verbose)
        else:
           MRtrixToSbiaFodRun(data_im,outFile,verbose)
                     
    except Exception,e:
        raise e

    if not os.path.isfile(outFile):
        raise SbiaException("Output file not found: " + outFile)

    return 0

#---------------------------------------------    

def index(l,m):
    """Returns the index in the basis ordering (0,0),(2,-2),(2,-1),(2,0)..."""
    return ((l+l**2)/2 + m)
    
#-----------------------------------------    
def sign(m):
    """Evaluates (-1)^m"""
    if m%2 == 0:
       return 1
    else:
       return -1
    
#--------------------------------------------    
def getLM(index):
    """Returns the (l,m) tuple from the index in the basis ordering (0,0),(2,-2),(2,-1),(2,0),..."""
    l=int((-1+sqrt(1+8*index))/2)
    m=index(l,0)-index
    return (l,m)
    
#---------------------------------------------------    

def SbiaFodToMRtrixRun(image, outFile, verbose=False):

    #mrtrix expects image dimensions of (x,y,z,a)
    d=image.getDataArray()[:,0,:,:,:]
    im=numpy.zeros(d.shape,'float32')
    
    #determine lmax
    lmax=( (8*len(im)+1)**.5 - 3 )/2
    
    if int(lmax) != lmax:
        raise SbiaException("Data file "+dataFile+" is not of a valid length!")
    lmax = int(lmax)

    if verbose:
        print "Performing change of basis..."
    for l in range(0,lmax+1,2): #gives [0,2,4,6,8]
        for m in range(-l,l+1):
            if m<0:
                im[index(l,m),:,:,:] = sqrt(2)/2 * d[index(l,-m),:,:,:]
            elif m==0:
                im[index(l,0),:,:,:] = d[index(l,0),:,:,:]
            elif m>0:
                im[index(l,m),:,:,:] = sign(m) * sqrt(2)/2 * d[index(l,-m),:,:,:]
            
    if verbose:
        print "Done!"
        
    #Nifti header stuff: we changed the axes
    im_header=image.asDict()
    extent=image.getExtent()
    im_header['dim']=[4,extent[0],extent[1],extent[2],extent[4],1,1,1]
     
    nim=NiftiImage(im,header=im_header)
    nim.save(outFile)
    if verbose:
        print "Output file created:  ",outFile
	
#----------------------------------------------------------	

def MRtrixToSbiaFodRun(image, outFile, verbose=False):

    #mrtrix expects image dimensions of (x,y,z,a)
    d=image.getDataArray()
    im=numpy.zeros(d.shape,'float32')
    
    #determine lmax
    lmax=( (8*len(im)+1)**.5 - 3 )/2
    
    if int(lmax) != lmax:
        raise SbiaException("Data file "+dataFile+" is not of a valid length!")
    lmax = int(lmax)

    if verbose:
        print "Performing change of basis..."
    for l in range(0,lmax+1,2): #gives [0,2,4,6,8]
        for m in range(-l,l+1):
            if m<0:
                im[index(l,m),0,:,:,:] = sign(abs(m)) * sqrt(2) * d[index(l,-m),:,:,:]
            elif m==0:
                im[index(l,0),0,:,:,:] = d[index(l,0),:,:,:]
            elif m>0:
                im[index(l,m),0,:,:,:] = sqrt(2) * d[index(l,-m),:,:,:]
            
    if verbose:
        print "Done!"
        
    #Nifti header stuff: we changed the axes
    im_header=image.asDict()
    extent=image.getExtent()
    im_header['dim']=[5,extent[0],1,extent[1],extent[2],extent[3],1,1]
    im_header['intent_code'] = 1007  

    nim=NiftiImage(im,header=im_header)
    nim.save(outFile)
    if verbose:
        print "Output file created:  ",outFile

main()
