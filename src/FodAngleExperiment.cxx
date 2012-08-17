#define VERBOSE_CXX
#define INFO_FILE

#include "itkSymRealSphericalHarmonicRep.h"
#include "itkGeodesicConcentrationFilter.h"
#include "itkGeneralFractionalAnisotropyImageFilter.h"
#include "itkVector.h"
#include "string.h"

#include "math.h"
#include "itkNiftiImageIO.h"
#include "spherePoints.h"
#include "FodAngleExperiment.h"
#include "iostream.h"
#include "itkNumericTraits.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_matrix.h"

/*
 * ComputeGCatSpherePoints.cxx
 *  
 * This executable computes Geodesic Concentration of RSH functions at 100 equally spaced
 * sphere points defined in "sphere1000points.h"
 *
 */



const unsigned int    NOrder = 8;
const unsigned int    Dimension = 3;
typedef float         ComponentType;

typedef itk::SymRealSphericalHarmonicRep< ComponentType, NOrder >                                       RshPixelType;
typedef ComponentType                                                                                   OutputPixelType;
typedef itk::Vector< ComponentType, Dimension >                                                         PeakPixelType;
typedef itk::Functor::GeodesicConcentration< RshPixelType, PeakPixelType, OutputPixelType >             GcCalcType;
typedef itk::Functor::GeneralFractionalAnisotropy< RshPixelType, OutputPixelType >                      GfaCalcType;

int main(int argc, char* argv[])
{
const int numVxls = 10;
GcCalcType gcCalc; 
GfaCalcType gfaCalc;

PeakPixelType first_peak, second_peak;
GcCalcType::OutputPixelType outputFirstGC[numVxls];
GcCalcType::OutputPixelType outputSecondGC[numVxls];
GfaCalcType::OutputPixelType outputGfa[numVxls];

RshPixelType   shCoefs[ numVxls ];
for ( int i = 0; i < numVxls; i ++ ) shCoefs[i] = crossesExp[i];


for( int i=0; i < numVxls; ++i ) 
{ 
  first_peak.SetElement(0, firstPeak[i][0]); 
  first_peak.SetElement(1, firstPeak[i][1]);
  first_peak.SetElement(2, firstPeak[i][2]);
  
  second_peak.SetElement(0, secondPeak[i][0]);
  second_peak.SetElement(1, secondPeak[i][1]);
  second_peak.SetElement(2, secondPeak[i][2]);
  
  //std::cout << dir << "\n";
  outputFirstGC[i] = gcCalc(shCoefs[i], first_peak);
  outputSecondGC[i] = gcCalc(shCoefs[i], second_peak);
  outputGfa[i] = gfaCalc(shCoefs[i]);
  //std::cout << outputFirstGC[i] << "\n";  
}

    #ifdef INFO_FILE
        FILE *fd = fopen("/sbia/home/raoarvin/FodAngleExpResults.txt", "w");
        //fprintf( fd, "\n\n\nRotated Coefficients Are:   \n\n");    
        
        for(int i=0; i< numVxls; ++i)
        {
           std::cout << outputFirstGC[i] << "\n";
           
           fprintf( fd, "%.12g\t" , outputFirstGC[i]);
           fprintf( fd, "%.12g\t" , outputSecondGC[i]);
           fprintf( fd, "%.12g\n" , outputGfa[i]);
        }
   
        fclose(fd);
    #endif
     
return EXIT_SUCCESS;
}//end main
