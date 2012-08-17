#define VERBOSE_CXX

#include "itkSymRealSphericalHarmonicRep.h"
#include "itkRshDerivativeFunctions.h"
#include "RshDerivativesTestData.h"
#include "string.h"
#include "iostream.h"
#include "itkReplaceSpecialFunctions.h"
#include "math.h"

typedef float                                                         ComponentType;

unsigned int getIndex(int l, int m)
{
  int b[5] = {0, 1, 6, 15, 28};
  return l+m+b[ l/2 ];
}

int main(int argc, char* argv[])
{


const unsigned int                                                             numCoeff = 45;
const unsigned int                                                             numTest   =  4;
const unsigned int                                                             Order     =  8;
const  ComponentType                                                  epsilon   = .0001; 
typedef itk::RshDerivativeFunctions<ComponentType, Order>             FunctionsClass; 
ComponentType                                                         temp;

ComponentType xValues[numTest] = {.25, .50, .75, .95 };
ComponentType del_x_PlmComputed[numTest][numCoeff];
ComponentType del_x2_PlmComputed[numTest][numCoeff];
FunctionsClass derivatives;

for( unsigned int k = 0; k < numTest; k++)
{
  for ( int l = 0; l <= int(Order); l+=2)
   {
     for( int m = -l; m <= l; m++)
     {
       //std::cout<< " P^l_m for " << "l=" << l << " and " << "m=" << m << " is: " << derivatives.del_x_LegendreP(l,m, 0.25) << "\n";
       // fill del_x and del_x2 arrays
       del_x_PlmComputed[k][ getIndex(l,m) ]= derivatives.del_x_LegendreP(l,m, xValues[k]);
       del_x2_PlmComputed[k][ getIndex(l,m) ]= derivatives.del_x2_LegendreP(l,m, xValues[k]);
     }
   }
}

for( unsigned int k = 0; k < numTest; k++)
{
  for ( int l = 0; l <= int(Order); l+=2)
   {
     for( int m = -l; m <= l; m++)
     { 
       temp = fabs( del_x_PlmComputed[k][ getIndex(l,m) ] -  del_x_PlmTest[k][ getIndex(l,m)]);
       if (  temp > epsilon )
           std::cout<< "Error!! Difference for dP^l_m at " << "l=" << l << " and " << "m=" << m << " is: " <<  temp << ", X value is " << xValues[k] << "\n";
           
       temp = fabs( del_x2_PlmComputed[k][ getIndex(l,m) ] -  del_x2_PlmTest[k][ getIndex(l,m)] );
       if (  temp > epsilon )
           std::cout<< "Error!! Difference for ddP^l_m at " << "l=" << l << " and " << "m=" << m << " is: " <<  temp << ", X value is " << xValues[k] << "\n";
 
     }
   }
}
return EXIT_SUCCESS;
}
