These files are provided and described for continuation of my work on
Geodesic Concentration and to facilitate further development of HARDI
clinical scalar measures. If you use them as suggested they should work
fine. Please read and understand the file you intend you use, to
determine if it suits your needs. For more information, please contact
me at `Arvind Rao`_

MAPLE 14 Files 
--------------

MAPLE 14 was used for various tasks: visualization of RSH functions,
Solid Angle-ODF and FOD modeling fitting, experiments.

- **EquidistantSpherePoints.mw** - Computes the 


Streamline Related Scripts 
--------------------------

*Some of these scripts may require the `SPLINE FIT`_ library*

-  **arclengthPP.m** - Computes the arclength of a curve represented as
   piecewise polynomial, given a time step and 1-d vector of sample time
   points, from the domain of the curve.

-  **CutStreamLine.m** - Given a streamline, represented by a series of
   3-d points, this script cuts the streamline based on cutting planes
   determined anatomically. Also resamples the streamline into a fixed
   number of points. This script was used specifically to cut
   streamlines so that all streamlines start and stop along common
   cross-sections of a particular fiber.

-  **truncateStreamLines.m** - Loop through streamlines output from
   MrTrix, and cuts streamlines (calls CutStreamLine.m). Each processed
   streamline is output into a csv format array of 3-d points, readable
   by ComputeScalarsOnStreamLine.cxx (See HardiScalars project).

-  **PlotStreamLinesWithColors.m** - Given a directory with truncated
   streamline files, mean and variance across corresponding points are
   computed and plotted for GC & GFA.
   
-  **write_mrtrix_track.m** - Writes streamline files compatible with
   MrTrix (i.e. tracks.tck). [written by unknown]
   
-  **read_mrtrix_track.m** - Reads streamline files (*.tck) produced by 
   MrTrix. [written by unknown]

For Statistical Analysis 
------------------------

**MaxClusterSizeStat.m** - Computes the size of the largest cluster in a
Hotelling T^2 statistical map, generated from permuted labels. This
function is dependent on the `SurfStat`_ library. This function can
easily be modified to compute the max cluster size of an univariate
t-statistic. The linear model considered is Y = 1 + Group + Age.

**MaxClusterSizeStat\_core.m** - Computes the max cluster size in the
same way as MaxClusterSizeStat.m, however it can be mcc’ed. It takes an
integer, ‘count’, to label the output file, output directory, and an
input .mat file with variables containing the vector GC maps and
associated labels.

.. _Arvind Rao: mailto:arvind.sbia@gmail.com
.. _SPLINE FIT: http://www.mathworks.com/matlabcentral/fileexchange/13812-splinefit
.. _SurfStat: http://www.math.mcgill.ca/keith/surfstat/
