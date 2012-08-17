#define VERBOSE_CXX

#include "itkSymRealSphericalHarmonicRep.h"
#include "itkGeneralFractionalAnisotropyImageFilter.h"
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
template <typename TComponentType, unsigned int NOrder, unsigned int ImageDim >
int runOverRshDim(std::string  rsh_file, std::string  output_file, bool normArg )
{

	typedef itk::Image< itk::SymRealSphericalHarmonicRep< TComponentType, NOrder > , ImageDim >             RshImageType; 
	typedef itk::Image< TComponentType, ImageDim >                                                          OutputImageType;
	typedef itk::GeneralFractionalAnisotropyImageFilter< RshImageType, OutputImageType >                    FilterType;
		                                                                                                
	typedef itk::ImageFileReader< RshImageType >                                                            RshReaderType;
	typedef itk::ImageFileWriter< OutputImageType >                                                         WriterType;

	typename RshReaderType::Pointer rsh_reader = RshReaderType::New();

	typename WriterType::Pointer writer = WriterType::New();
	typename FilterType::Pointer filter = FilterType::New();

	if(normArg) 
	  filter->NormalizationOn(); //turn ON normalization. Filter default is OFF.

	rsh_reader->SetFileName( rsh_file ); 
	writer->SetFileName( output_file + "_gfa.nii.gz");

	try
	{
		std::cout<< "Calculating General Fractional Anisotropy Maps" <<"\n";
		filter->SetInput(rsh_reader->GetOutput());
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
		std::cout<< "Writing GFA Maps Image" <<"\n";
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

	try
	{
		CmdLine cmd("Calculate General Fractional Anisotropy given an RSH image", "HardiScalars", "This executable computes GFA","GeneralFractionalAnisotropy -a ${InputDir}/${ID}_FOD_Warp.nii.gz -p ${OutputDir}/${ID} -N 1", "1.0");

		StringArg rshfileArg("a","rshfile","RSH Image File. Expects a RSH images of pixel type itk::SymRealSphericalHarmonicRep and order ",true,"homer","string");
		cmd.add( rshfileArg );

		StringArg prefixArg("p","prefix","Prefix for Output Image. Output is named 'prefix_gfa.nii.gz'",true,"homer","string");
		cmd.add( prefixArg );

		ValueArg<bool> normalizationArg("N","normalize","Normalization of RSH pixeltype, done within the operator. Default value is FALSE. Input should be either 0/1. ",false,false,"bool");
		cmd.add( normalizationArg );

		//TCLAP::ValueArg<std::string> outputdirArg("o","outputdir","Output Directory. Defaults to current directory.",true,"homer","string");
		//cmd.add( outputdirArg );

		// Parse the argv array.
		cmd.parse( argc, argv );

		readerRsh->SetFileName( rshfileArg.getValue() ); 
	try 
	 { 
	   readerRsh->Update();
	 } 
	catch( itk::ExceptionObject & err ) 
	 { 
	   std::cerr << "ExceptionObject caught !" << std::endl; 
	   std::cerr << err << std::endl; 
	   return EXIT_FAILURE;
	 } 
	      
	unsigned int RshDim=readerRsh->GetOutput()->GetVectorLength();

	#ifdef VERBOSE_CXX 
	std::cout<< "RSH File Name: " << rshfileArg.getValue() << "\n";
	std::cout<< "Prefix For Filename: " << prefixArg.getValue() << "\n";
	//std::cout<< "Output Directory: " << outputdirArg.getValue() << "\n";
	std::cout<< "length of RSH vector images: "<< RshDim <<"\n"; 
	std::cout<< "value of Normalization Arg: " << normalizationArg.getValue() << "\n";
	#endif

	switch( RshDim )
	{
		case 6:
		    return runOverRshDim< ComponentType, 2, Dimension >(rshfileArg.getValue(), prefixArg.getValue(), normalizationArg.getValue() );
		    break;
		case 15:
		    return runOverRshDim< ComponentType, 4, Dimension >(rshfileArg.getValue(),  prefixArg.getValue(), normalizationArg.getValue() );
		    break;
		case 28:
		    return runOverRshDim< ComponentType, 6, Dimension >(rshfileArg.getValue(), prefixArg.getValue(), normalizationArg.getValue() );
		    break;
		case 45:
		    return runOverRshDim< ComponentType, 8, Dimension >(rshfileArg.getValue(), prefixArg.getValue(), normalizationArg.getValue() );
		    break;    
		default:
		    return runOverRshDim< ComponentType, 8, Dimension >(rshfileArg.getValue(), prefixArg.getValue(), normalizationArg.getValue() );
	}

	}
	catch(ArgException &e)  // catch any exceptions
	{ 
	    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
	}

       return EXIT_SUCCESS;
}//end main
