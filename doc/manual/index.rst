HardiScalars
============

Introduction 
------------

HARDI models, such as FOD and ODF, provide added information about water diffusion in areas of complex white matter (CWM), i.e 
fiber crossings. To take advantage of this extra information and facilate standard voxel-based analysis, new measures are
needed. HardiScalars primarily implements a directed diffusion measure called Geodesic Concentration (GC) and General Fractional
Anisotropy (GFA). GC is defined as the concentration of a HARDI model in a specific direction. As such, there are as many 
GC values as there are fiber populations represented by the HARDI model. Further details about this method are contained 
in our `MMBIA paper`_. 

The majority of the code contained in this project is centered around experiments written during the course of investigating 
properties of GC. There are various image analysis tools like *ComputeAngleBtwPeaks.py* which aligns subject primary and secondary
peaks with those of the template, and scripts for producing statistical output like *gc_norm.sh*. Statistical and processing scripts 
are included as examples.


Running Environment
-------------------

- GC and GFA are command-line software tools that run on Linux (Ubuntu, Redhat, CentOS) and Mac OS.
- Image proccessing scripts and scripts for producing statistical output are also command-line tools. 
- `MAPLE and MATLAB files`_ require MAPLE and MATLAB environments respectively.    


Features
--------

- Implementations of clinical scalar measures for HARDI: Geodesic Concentration, General Fractional Anisotropy.

- /src/statistic folder contains BASH scripts which implement group studies conducted for GC and GFA.  

- MATLAB files which interpolate, smooth, and truncate MrTrix streamlines.

- MAPLE files for visualization of ODF/FOD, implementations of Solid Angle ODF and CSD modeling of DW signals 

- C++ code for smoothing of FOD images (written by Harini Eavani)

- C++ code for computation of Solid Angle images from DW data. 


Applications
------------

- Voxel-based Analysis of HARDI images.

- C++ executable for computing GC and GFA along streamlines. In combination with MATLAB scripts provided, it 
  can be used as a basis for tract based statistics. 


.. _Arvind Rao: mailto:arvind.sbia@gmail.com
.. _MMBIA paper: http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6164756
.. _MAPLE and MATLAB files: matlab_maple_scripts.rst
