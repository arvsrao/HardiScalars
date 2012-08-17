#include "itkSymRealSphericalHarmonicRep.h"
#include "itkPeakFindingCalculator.h"

/*****
*
* PeakPixelCalculator.cxx
*
* Uses the Peak-Finder to compute peaks for one
* RSH pixel. Outputs all 9 local maxima found.:w
*
* User must provide the RSH coefficients of the pixel. 
*
*
*/


/****************************************************
                       MAIN
****************************************************/

int main(int argc, char* argv[])
{

const unsigned int    NOrder = 8; 
const unsigned int    Dimension = 3;

typedef float         ComponentType;

typedef itk::SymRealSphericalHarmonicRep< ComponentType, NOrder>                        PixelType;
typedef itk::PeakFindingCalculator< PixelType, ComponentType, ComponentType >           PeakFindingCalculator;
typedef PeakFindingCalculator::InputType                                                PeakFindingInputType; 

PeakFindingCalculator::Pointer pFinder = PeakFindingCalculator::New(); 
std::vector< ComponentType >                                                            odPeaks;
PeakFindingInputType                                                                    C; 
 
//ODF image is read

if(argc < 44)
 {
    std::cerr << " Usage:\n" << "\n";
    std::cerr << "Coefficients of RSH function (must be length 45) " <<"\n";
    return 0; 
 }
else 
 {  
    std::cout<< " Enter the RSH coefficients of your function (must length 45):" << "\n\n";
    for(int i=0; i<45; i++) C[i] = float( std::atof(argv[i]) );  
 }

odPeaks.clear();
pFinder->SetCoefficients( C );
pFinder->SetPeakFiltering( 1 );
pFinder->SetAcceleration( 0 );
pFinder->GetRawPeaks( odPeaks );

for(int i=0; i < 36; i++)
 {
    if( i % 4 ==0 )
     {
       std::cout<< "Peak " << float(i)/float(4) + 1 <<" :   ";
       std::cout<< "("<< odPeaks[i] << ", " << odPeaks[i+1] << ", " << odPeaks[i+2] << ")\n";
     }  
 }
 
//for(int i=0; i<45; i++) std::cout <<" "<< C[i]; 

for(int i=0; i<36; i++) std::cout <<" "<< odPeaks[i];

return 0;
}//end main
