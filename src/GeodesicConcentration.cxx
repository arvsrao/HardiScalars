#define VERBOSE_CXX

#include "itkSymRealSphericalHarmonicRep.h"
#include "itkGeodesicConcentrationFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVector.h"
#include "string.h"
#include "itkNiftiImageIO.h"
#include "iostream.h"
#include "basis.h"

using namespace sbia::basis;

//Template parameters are only handled at compile time. 
template <typename TComponentType, unsigned int NOrder, unsigned int TPeakDimension, unsigned int ImageDim >
int runOverRshDim(std::string  rsh_file,  std::string  peak_file,  std::string  output_file, bool normArg )
{

typedef itk::Image< itk::SymRealSphericalHarmonicRep< TComponentType, NOrder > , ImageDim >             RshImageType; 
typedef itk::Image< itk::Vector< TComponentType, TPeakDimension >, ImageDim >                           VecImageType;
typedef itk::Image< TComponentType, ImageDim >                                                          OutputImageType;
typedef itk::GeodesicConcentrationFilter< RshImageType, VecImageType, OutputImageType >                 FilterType;
                                                                                                        
typedef itk::ImageFileReader< VecImageType >                                                            PeaksReaderType;
typedef itk::ImageFileReader< RshImageType >                                                            RshReaderType;
typedef itk::ImageFileWriter< OutputImageType >                                                         WriterType;

typename RshReaderType::Pointer rsh_reader = RshReaderType::New();
typename PeaksReaderType::Pointer peaks_reader = PeaksReaderType::New();

typename WriterType::Pointer writer = WriterType::New();
typename FilterType::Pointer filter = FilterType::New();

if(normArg) 
   filter->NormalizationOn(); //turn on normalization.

rsh_reader->SetFileName( rsh_file ); 
peaks_reader->SetFileName( peak_file );
writer->SetFileName( output_file + "_gc.nii.gz");

try
{
	std::cout<< "Calculating Geodesic Concentration Maps" <<"\n";
	filter->SetInput1(rsh_reader->GetOutput());
	filter->SetInput2(peaks_reader->GetOutput());
	filter->Update();
}
catch( itk::ExceptionObject & err )
{
	std::cerr <<"ExceptionObject caught in the Filter !" <<"\n";
	std::cerr << err << std::endl;
	return EXIT_FAILURE;
}

try
{
	std::cout<< "Writing GC Map Image" <<"\n";
	writer->SetImageIO( itk::NiftiImageIO::New() );
	writer->SetInput( filter->GetOutput() );
	writer->Update();
}
catch( itk::ExceptionObject & err )
{
	std::cerr << "ExceptionObject caught in the Writer !" << std::endl;
	std::cerr << err << std::endl;
	return EXIT_FAILURE;
}

return EXIT_SUCCESS;
}//end runOverRshDim.


int main(int argc, char* argv[])
{

const unsigned int    Dimension = 3;
typedef float         ComponentType;

typedef itk::VectorImage< ComponentType, Dimension >        VectorImageType;
typedef itk::ImageFileReader< VectorImageType >             ReaderType;
ReaderType::Pointer readerRsh = ReaderType::New();
ReaderType::Pointer readerPeaks = ReaderType::New();

try
{
        CmdLine cmd("Calculate Geodesic Concentration given an image of peaks and an RSH image", "HardiScalars", "This executable computes Geodesic Concentration",
        "GeodesicConcentration -N 1 -a ${RshDir}/${ID}_FOD_Warp_Norm.nii.gz -b ${OutputDir}/${ID}_Norm_FOD_Warp_relabeled_1st_peak.nii.gz -p ${OutputDir}/${ID}_relabeled_1st_peak", "1.0");


	StringArg rshfileArg("a","rshfile","RSH Image File. Expects a RSH images of pixel type itk::SymRealSphericalHarmonicRep and order ",true,"homer","string");
	cmd.add( rshfileArg );

	StringArg peakfileArg("b","peakfile","Peak Image File. Image of 3-component vectors, representing directions in 3-d real euclidean space.",true,"homer","string");
	cmd.add( peakfileArg );

	StringArg prefixArg("p","prefix","Prefix for Output Image. Output is named 'prefix_gc.nii.gz'",true,"homer","string");
	cmd.add( prefixArg );

	ValueArg<bool> normalizationArg("N","normalize","Normalization of RSH pixeltype, done within the operator. Default value is FALSE. Input should be either 0/1. ",false,false,"bool");
	cmd.add( normalizationArg );

	//TCLAP::ValueArg<std::string> outputdirArg("o","outputdir","Output Directory. Defaults to current directory.",true,"homer","string");
	//cmd.add( outputdirArg );

	// Parse the argv array.
	cmd.parse( argc, argv );

	readerRsh->SetFileName( rshfileArg.getValue() ); 
	readerPeaks->SetFileName( peakfileArg.getValue() );
	try 
	 { 
	   readerRsh->Update();
	   readerPeaks->Update(); 
	 } 
	catch( itk::ExceptionObject & err ) 
	 { 
	   std::cerr << "ExceptionObject caught !" << std::endl; 
	   std::cerr << err << std::endl; 
	   return EXIT_FAILURE;
	 } 
	    
	unsigned int PeaksDim=readerPeaks->GetOutput()->GetVectorLength();    
	unsigned int RshDim=readerRsh->GetOutput()->GetVectorLength();

	#ifdef VERBOSE_CXX 
	std::cout<< "RSH File Name: " << rshfileArg.getValue() << "\n";
	std::cout<< "Peak Image Filename: " << peakfileArg.getValue()  << "\n";
	std::cout<< "Prefix For Filename: " << prefixArg.getValue() << "\n";
	//std::cout<< "Output Directory: " << outputdirArg.getValue() << "\n";
	std::cout<< "length of peak vectors: " << PeaksDim <<"\n";
	std::cout<< "length of RSH vector images: "<< RshDim <<"\n"; 
	#endif

	if( PeaksDim != 3 )
	 {
	    std::cerr<< " only 1 direction is supported at this time!!\n " << "\n";
	    return EXIT_FAILURE;
	 }

	switch( RshDim )
	{
	case 6:
	    return runOverRshDim< ComponentType, 2, 3, Dimension >(rshfileArg.getValue(), peakfileArg.getValue(), prefixArg.getValue(), normalizationArg.getValue() );
	    break;
	case 15:
	    return runOverRshDim< ComponentType, 4, 3, Dimension >(rshfileArg.getValue(), peakfileArg.getValue(), prefixArg.getValue(),normalizationArg.getValue()  );
	    break;
	case 28:
	    return runOverRshDim< ComponentType, 6, 3, Dimension >(rshfileArg.getValue(), peakfileArg.getValue(), prefixArg.getValue(), normalizationArg.getValue() );
	    break;
	case 45:
	    return runOverRshDim< ComponentType, 8, 3, Dimension >(rshfileArg.getValue(), peakfileArg.getValue(), prefixArg.getValue(),normalizationArg.getValue()  );
	    break;    
	default:
	    return runOverRshDim< ComponentType, 8, 3, Dimension >(rshfileArg.getValue(), peakfileArg.getValue(), prefixArg.getValue(), normalizationArg.getValue() );
}

}
catch(ArgException &e)  // catch any exceptions
{ 
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
}

return EXIT_SUCCESS;
}//end main
