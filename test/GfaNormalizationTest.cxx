#define VERBOSE_CXX
#define INFO_FILE

#include "itkSymRealSphericalHarmonicRep.h"
#include "itkGeodesicConcentrationFilter.h"
#include "itkGeneralFractionalAnisotropyImageFilter.h"

#include "math.h"
#include "vector.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVector.h"
#include "string.h"
#include "itkNiftiImageIO.h"
#include "tclap/CmdLine.h"
#include "iostream.h"


//print STL vector
template< typename TPixelType>
void PrintVector( std::vector<TPixelType> &vec )
{
   for(int j=0; j < int(vec.size()); j++) std::cout << vec[j] << "\n";
}


//Template parameters are only handled at compile time. 
template <typename TComponentType, unsigned int NOrder, unsigned int ImageDim >
int runOverRshDim(std::string  rsh_file, std::vector<int> vec )
{

typedef itk::SymRealSphericalHarmonicRep< TComponentType, NOrder >                                      RshPixelType;
typedef itk::Image< RshPixelType, ImageDim >                                                            RshImageType;   
typedef itk::ImageFileReader< RshImageType >                                                            RshReaderType;

typename RshReaderType::Pointer rsh_reader = RshReaderType::New();
typedef typename RshImageType::IndexType                                                                RshIndexType;
typedef itk::Functor::GeneralFractionalAnisotropy< RshPixelType, TComponentType >                       GfaCalcType;

//std::cout<< "Voxel location is: " << vec << "\n";
try
{
rsh_reader->SetFileName( rsh_file ); 
//peaks_reader->SetFileName( peak_file );
rsh_reader->Update();
//peaks_reader->Update();
}
catch( itk::ExceptionObject & err )
{
std::cerr <<"ExceptionObject caught in the Filter !" <<"\n";
std::cerr << err << std::endl;
return EXIT_FAILURE;
}

//Actual experiment...
GfaCalcType gfaCalc; 
RshIndexType  voxel;

for(unsigned int j=0; j < ImageDim; j++) voxel[j] = vec[j];
std::cout<< "Voxel location is: " << voxel << "\n";


RshImageType *rshImg = rsh_reader->GetOutput();
RshPixelType rsh = rshImg->GetPixel(voxel);

std::cout<< "coefficients of Un-normalized FOD" << rsh << "\n";
std::cout<< "Un-normalized GFA is: " << gfaCalc( rsh ) << "\n";

rsh.Normalize();
std::cout<< "coefficients of Normalized FOD" << rsh << "\n";
std::cout<< "Normalized GFA is: " << gfaCalc( rsh ) << std::endl;

return EXIT_SUCCESS;
}//end runOverRshDim.

int main(int argc, char* argv[])
{

const unsigned int    Dimension = 3;
typedef float         ComponentType;

typedef itk::VectorImage< ComponentType, Dimension >                                                  VectorImageType;
typedef itk::ImageFileReader< VectorImageType >                                                       ReaderType;
typedef itk::Image< itk::Vector<ComponentType, Dimension>, Dimension>                                 PeakImageType;

ReaderType::Pointer readerRsh = ReaderType::New();

try
{

TCLAP::CmdLine cmd("Calculate Geodesic Concentration & GFA for a streamline", ' ', "1.0");

TCLAP::ValueArg<std::string> rshfileArg("d","rshfile","RSH Image File. Expects a RSH images of pixel type itk::SymRealSphericalHarmonicRep and order ",true,"homer","string");
cmd.add( rshfileArg );

TCLAP::MultiArg<int> pointArg("v", "coordinate", "x,y,z voxel coordinates", true, "int" );
cmd.add( pointArg );

// Parse the argv array.
cmd.parse( argc, argv );

try 
 { 
   readerRsh->SetFileName( rshfileArg.getValue() ); 
   readerRsh->Update();
  
  //readerPeaks->SetFileName( peakfileArg.getValue() );
  // readerPeaks->Update(); 
 } 
catch( itk::ExceptionObject & err ) 
 { 
   std::cerr << "ExceptionObject caught !" << std::endl; 
   std::cerr << err << std::endl; 
   return EXIT_FAILURE;
 }   

const unsigned int RshDim=readerRsh->GetOutput()->GetVectorLength();

#ifdef VERBOSE_CXX 
std::cout<< "RSH File Name: " << rshfileArg.getValue() << "\n";

std::cout<< "length of RSH vector images: "<< RshDim <<"\n"; 
#endif

switch( RshDim )
{
case 6:
    return runOverRshDim< ComponentType, 2, Dimension >(rshfileArg.getValue(), pointArg.getValue() );
    break;
case 15:
    return runOverRshDim< ComponentType, 4, Dimension >(rshfileArg.getValue(), pointArg.getValue() );
    break;
case 28:
    return runOverRshDim< ComponentType, 6, Dimension >(rshfileArg.getValue(), pointArg.getValue() );
    break;
case 45:
    return runOverRshDim< ComponentType, 8, Dimension >(rshfileArg.getValue(), pointArg.getValue() );
    break;    
default:
    return runOverRshDim< ComponentType, 8, Dimension >(rshfileArg.getValue(), pointArg.getValue() );
}
                                                                                     
}
catch(TCLAP::ArgException &e)  // catch any exceptions
{ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
}

return EXIT_SUCCESS; 
}//end main
