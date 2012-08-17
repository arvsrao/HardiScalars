##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of software.
#
# Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
# See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <sbia-software at uphs.upenn.edu>
##############################################################################

# ============================================================================
# library target(s)
# ============================================================================

# Add library target for each library using basis_add_library().
#
# This command can not only be used to build libraries from C/C++, but also source
# code written in other programming languages such as Java, Python, Perl, BASH,
# and MATLAB. Note that here we consider modules written in a scripting language
# which are no executables but to be included by other scripts written in the
# particular language as libraries.



# ============================================================================
# executable target(s)
# ============================================================================

# Add executable target for each executable program using basis_add_executable().
#
# This command can not only be used to build executables from C/C++, but also
# source code written in other programming languages such as Java, Python, Perl,
# BASH, and MATLAB.

# ----------------------------------------------------------------------------
# add executable which links to the ITK and the SBIA Toolkit
macro (hardiscalars_add_executable SRC)
  basis_get_source_target_name(TGT "${SRC}" NAME_WE)
  basis_add_executable(${TGT} LIBEXEC ${SRC} ${ARGN})
  basis_target_link_libraries(${TGT} ITKCommon ITKIO ITKNumerics ${SBIATOOLKIT_LIBS})
endmacro ()

# ----------------------------------------------------------------------------
# list of executables implemented in C++
set (EXES

)

# ----------------------------------------------------------------------------
# list of executables implemented in scripting language
set (SCRIPTS

  compute2ndPeakMask.sh
  ConvertToMrtrixBasis.py
  MakeODFs.sh
  MakeSolidAngleODFs.sh
  PeakMapsMakingSA.sh
  PeakMapsMaking.sh
  register_wm_masks.sh
  registration.sh
  SmoothingTest.sh
  SplitGcMaps.py
  SplitVectorImage.py  
)

# ----------------------------------------------------------------------------
# add executables to build configuration
foreach (EXE IN LISTS EXES)
  hardiscalars_add_executable(${EXE})
endforeach ()

foreach (SCRIPT IN LISTS SCRIPTS)
  basis_add_executable(${SCRIPT})
endforeach ()
