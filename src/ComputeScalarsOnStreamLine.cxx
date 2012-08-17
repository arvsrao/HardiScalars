//#define VERBOSE_CXX
//#define INFO_FILE
//#define PRINT_VXLS
//#define TESTING_CXX

#include "itkSymRealSphericalHarmonicRep.h"
#include "itkGeodesicConcentrationFilter.h"
#include "itkGeneralFractionalAnisotropyImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "itkContinuousIndex.h"
#include "itkPointSet.h"

#include "math.h"
#include "vector.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVector.h"
#include "string.h"
#include "itkNiftiImageIO.h"
#include "iostream.h"
#include "fstream.h"

#include "basis.h"


using namespace sbia::basis;


//print STL vector
template< typename TPixelType>
void PrintVector( std::vector<TPixelType> &vec )
{
   for(int j=0; j < int(vec.size()); j++) std::cout << vec[j] << "\n";
}

/*
*  Change from image space coordinates to voxel coordinates
*    
*
*  Taken from qoffest of the template subject, which in the case of this
*  experiment is R0090 (Tobacco Study). 
*/
template<typename TComponentType, typename TPixelType>
void ChangeCoordinates( std::vector<TPixelType> &points )
{
       TPixelType qoffset;
       qoffset[0] = -124.803101;
       qoffset[1] = -136.537231;
       qoffset[2] = -84.061363;
       
      for(int j=0; j < int(points.size()); j+=1) 
       {
	      points[j] = (points[j] - qoffset) / TComponentType(2);
       } 
}

//Read the Stream File.
template< typename TComponentType, typename TPixelType>
void ReadStreamFile( std::string filename, std::vector<TPixelType> &points )
{
   TPixelType dir; 
   float num;
   std::string line, temp;
   std::vector<TComponentType> vec;
    
   std::ifstream indata( filename.c_str() ); // opens the file
   
   
   if(!indata) //test to see if the file has been openned successfully
   { 
      cerr << "Error: file could not be opened" << endl;
      exit(1);
   }
   else 
   {
       	   while( getline(indata, line) )
	    {
        
	         std::string::size_type first = line.find(',');
	         num = atof( line.substr(0, first).c_str() ); 
	         vec.push_back( num );
	     
	     #ifdef DEBUG 
	         std::cout << "first comma found at: " << int(first) << "\n";
	     #endif  
	           
	         std::string::size_type second = line.substr(first+1, line.size()-1).find(',');
	         num = atof( line.substr(first+1, first+second+1).c_str() );
	         vec.push_back( num );
	    
	      #ifdef DEBUG  
	         std::cout << "second comma found at: " << int(first + 1 + second) << "\n";
	      #endif
	         
	         num = atof( line.substr(first+second+2, line.size()-2).c_str() );                
	         vec.push_back( num );
	     
	      #ifdef DEBUG
	         std::cout << "last comma found at: " << int(line.size()-2) << "\n";
	      #endif
	         
	    }
	   indata.close();   
	   
	   #ifdef TESTING_CXX
	      std::cout<< "Printing the vector internal to ReadStreamFile(...)" << "\n";
	      PrintVector<TComponentType>( vec );
	   #endif
	    
	   for(int j=0; j < int(vec.size()); j+=3) 
	   {
	      dir[0] = vec[j];
	      dir[1] = vec[j+1];
	      dir[2] = vec[j+2];
	      points.push_back( dir );
	   }
    }	   
}


template< typename TComponentType >
void WriteScalarValuesFile( std::vector<TComponentType> &parm, std::vector<TComponentType> &gfaVals, std::vector<TComponentType> &gcVals, std::string filename )
{
           std::ofstream outdata( filename.c_str() ); // opens the file
         
          if( !outdata ) //test to see if the file has been openned successfully
	   { 
	      std::cerr << "Error: file could not be opened" << endl;
	      exit(1);
	   }
	  else 
	   {
	       //   outdata << "Time \t GFA \t GC \n";
	       for(unsigned int j= 0; j < parm.size(); ++j)
	       {
	         // outdata << parm[j] << " ";
	          outdata << gfaVals[j] << " ";
	          outdata << gcVals[j]  << "\n";
	       }
	       outdata.close();
	   }    
}


//Template parameters are only handled at compile time. 
template< typename TComponentType, typename TPixelType, unsigned int NOrder, unsigned int ImageDim >
int runOverRshDim(std::string  rsh_file, std::string out_filename, std::vector<TPixelType> &points)
{
const unsigned int ParametricDim = 1; 

typedef itk::SymRealSphericalHarmonicRep< TComponentType, NOrder >                                      RshPixelType;
typedef itk::Image< itk::Vector<TComponentType, ImageDim>, ImageDim>                                    PeakImageType;
typedef itk::Image< RshPixelType, ImageDim >                                                            RshImageType; 

typedef itk::ImageFileReader< RshImageType >                                                            RshReaderType;

typedef typename RshImageType::IndexType                                                                RshIndexType;
typedef typename PeakImageType::PixelType                                                               PeakPixelType;

typedef itk::Image< PeakPixelType, ParametricDim>                                                      CurveImageType;
typedef itk::PointSet< typename CurveImageType::PixelType, ParametricDim>                              PointSetType;
typedef typename PointSetType::PointType                                                               PointType; 

typedef itk::BSplineScatteredDataPointSetToImageFilter< PointSetType, CurveImageType>                  CurveInterpolatorType; 
typedef typename CurveInterpolatorType::GradientType                                                   GradientType;    

typedef itk::Functor::GeodesicConcentration< RshPixelType, PeakPixelType, TComponentType >             GcCalcType;
typedef itk::Functor::GeneralFractionalAnisotropy< RshPixelType, TComponentType >                      GfaCalcType;
typedef itk::VectorLinearInterpolateImageFunction< RshImageType , TComponentType >                     InterpolationType;
typedef typename InterpolationType::OutputType                                                         InterpolatorOutputType;

//set image readers
typename RshReaderType::Pointer rsh_reader = RshReaderType::New();
typename PointSetType::Pointer pointSet = PointSetType::New();
typename InterpolationType::Pointer interpolator = InterpolationType::New();
typename CurveInterpolatorType::Pointer curveInterpolator = CurveInterpolatorType::New();

try
{
    rsh_reader->SetFileName( rsh_file ); 
    rsh_reader->Update();
    interpolator->SetInputImage( rsh_reader->GetOutput() );
}
catch( itk::ExceptionObject & err )
{
    std::cerr <<"ExceptionObject caught in the Filter !" <<"\n";
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
}

// Declare Scalar calculators and variables...
//RshImageType *rshImg = rsh_reader->GetOutput();
std::vector<TPixelType> voxels = points;

//for(int j=0; j < int(points.size()); j++) voxels.push_back( points[j] ); 

ChangeCoordinates<TComponentType, TPixelType>( voxels );
//Print the vector.
#ifdef PRINT_VXLS
  PrintVector<TPixelType>( voxels );
#endif


// Instantiate GC calculator with normalization.
GcCalcType gcCalc;
gcCalc.NormalizationOn();

GfaCalcType gfaCalc; 
itk::Vector<TComponentType, 3> tempVec;
RshPixelType rsh; 

/*
*   Fill points vector;
*
*
*/

for(unsigned int i = 0; i < voxels.size(); ++i)
{
    PointType point;
    point[0] = TComponentType(i) / TComponentType( voxels.size() );
    //TPixelType tmpVxl;
    pointSet->SetPoint(i, point);
    pointSet->SetPointData(i, voxels[i]);
}

/*
* Compute tangents at points; they are peak directions. 
* Parametrize the curve however is convenient, because 
* the tangent obtained is normalized inside the GC filter, anyway. 
*
*/
typename CurveImageType::SpacingType spacing; spacing.Fill( 1.0 );  //voxel spacing
typename CurveImageType::SizeType size; size.Fill( voxels.size() );
typename CurveImageType::PointType origin; origin.Fill( 0.0 );

curveInterpolator->SetSize( size );
curveInterpolator->SetOrigin( origin );
curveInterpolator->SetSpacing( spacing);
curveInterpolator->SetInput( pointSet );

curveInterpolator->SetSplineOrder( 3 );
typename CurveInterpolatorType::ArrayType ncps;
ncps.Fill( 4 );

curveInterpolator->SetNumberOfControlPoints( ncps );
curveInterpolator->SetNumberOfLevels( 5 );    
curveInterpolator->SetGenerateOutputImage( false );

try
{
    curveInterpolator->Update();
}
catch(...) //ellipse handles catch any exception no matter what kind is thrown
{
    std::cerr <<"ExceptionObject caught in the Curve Interpolator !!" <<"\n";
    return EXIT_FAILURE;
}

std::vector<TComponentType> time, gfaVals, gcVals;

//Compute tangent lines at points. 
for( unsigned int i=0; i < voxels.size(); ++i)
{
    GradientType gradient;
    PointType ind;     
    InterpolatorOutputType vec;
    RshPixelType rsh; 
    itk::ContinuousIndex<TComponentType, 3> pixel;

    tempVec = voxels[i];
    time.push_back( TComponentType(i) / TComponentType( voxels.size() ) );
  
        
    for(unsigned int j  = 0; j < ImageDim; j++) 
    {
      ind[j] = tempVec[j];    
      pixel[j] = tempVec[j];
    }
    curveInterpolator->EvaluateGradientAtPoint( ind, gradient );
    vec = interpolator->EvaluateAtContinuousIndex( pixel );
 #ifdef DEBUG
    std::cout<< "# of rows " << gradient.Rows() <<"\n";
    std::cout<< "# of cols " << gradient.Cols() <<"\n";
 #endif   
    //Cast to Gradient to PeakPixelType
    PeakPixelType tempPeak;
    for(unsigned int j = 0; j < ImageDim; j++) tempPeak[j] = gradient(j, 0); 
   
    //Castmaek to RSH array type.
    for(int i = 0; i < RshPixelType::Dimension; ++i) { rsh[i] = vec[i]; }
    
    gfaVals.push_back( gfaCalc( rsh ) );
    gcVals.push_back( gcCalc(rsh, tempPeak ) ); 
   
#ifdef DEBUG
    std::cout<< "gradient at " << voxels[i] << " is: " << gradient << "\n";
#endif
}
    
WriteScalarValuesFile<TComponentType>( time, gfaVals, gcVals, out_filename );

#ifdef DEBUG
   itk::ContinuousIndex<TComponentType, 3> pixel;
   itk::Vector<TComponentType, 3> blah = voxels[10]; 
   for(int i  = 0; i < 3; ++i) pixel[i] = blah[i];
   interpolator->SetInputImage( rsh_reader->GetOutput() );
   RshPixelType something = interpolator->EvaluateAtContinuousIndex( pixel );
   for(int j = 0; j < 45; ++j) std::cout<< something[j] << "\n";
#endif

//for( int j=0; j < int(voxel.size()); j++)
//{
//    RshPixelType rsh = rshImg->GetPixel( voxels[j] );
//    
//}

//std::cout<< "GFA is: " << gfaCalc( rsh ) << "\n";
//std::cout<< "GC is: " << gcCalc( rsh , dir) << "\n";

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
        CmdLine cmd("Calculate Geodesic Concentration & GFA for a streamline", "HardiScalars", "This executable computes Geodesic Concentration & GFA for a streamline",
        "ComputeScalarsOnStreamLine -p $streamfile -d $template -o ./scalarOutput.txt","1.0");

        StringArg rshfileArg("d","rshfile","RSH Image File. Expects a RSH images of pixel type itk::SymRealSphericalHarmonicRep and order ",true,"homer","string");
        StringArg streamArg("p","streamfile","Stream-Line Image File, representing a fiber tract in the image coordinate system.",true,"homer","string");
        StringArg outputArg("o","outfile","Name of output file name",true,"homer","string");

        cmd.add( rshfileArg );
        cmd.add( streamArg );
        cmd.add( outputArg );

        //TCLAP::MultiArg<int> pointArg("v", "coordinate", "x,y,z voxel coordinates", true, "int" );
        //cmd.add( pointArg );

        // Parse the argv array.
        cmd.parse( argc, argv );

        try 
         { 
           readerRsh->SetFileName( rshfileArg.getValue() ); 
           readerRsh->Update();
         } 
        catch( itk::ExceptionObject & err ) 
         { 
           std::cerr << "ExceptionObject caught !" << std::endl; 
           std::cerr << err << std::endl; 
           return EXIT_FAILURE;
         }   

        const unsigned int RshDim=readerRsh->GetOutput()->GetVectorLength();

        //Read the input file
        std::vector<PeakImageType::PixelType> points, voxels; 
        try
        {
            ReadStreamFile<ComponentType, PeakImageType::PixelType>( streamArg.getValue(), points );
        }
        catch (std::bad_alloc& ba)
        {
             std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
            return 0;
        }

        //Print the vector.
        #ifdef VERBOSE_CXX
          PrintVector<PeakImageType::PixelType>( points );
        #endif

        #ifdef VERBOSE_CXX 
          std::cout<< "RSH File Name: " << rshfileArg.getValue() << "\n";
          std::cout<< "Stream-Line Filename: " << streamArg.getValue()  << "\n";
          std::cout<< "length of RSH vector images: "<< RshDim <<"\n"; 
        #endif

        switch( RshDim )
        {
            case 6:
                return runOverRshDim< ComponentType, PeakImageType::PixelType, 2, Dimension >(rshfileArg.getValue(), outputArg.getValue(), points );
                break;
            case 15:
                return runOverRshDim< ComponentType, PeakImageType::PixelType, 4, Dimension >(rshfileArg.getValue(), outputArg.getValue(), points );
                break;
            case 28:
                return runOverRshDim< ComponentType, PeakImageType::PixelType, 6, Dimension >(rshfileArg.getValue(), outputArg.getValue(), points );
                break;
            case 45:
                return runOverRshDim< ComponentType, PeakImageType::PixelType, 8, Dimension >(rshfileArg.getValue(), outputArg.getValue(), points );
                break;    
            default:
                return runOverRshDim< ComponentType, PeakImageType::PixelType, 8, Dimension >(rshfileArg.getValue(), outputArg.getValue(), points );
        }
                                                                                         
    }
    catch(ArgException &e)  // catch any exceptions
    { 
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
    }

    return EXIT_SUCCESS; 
}//end main
