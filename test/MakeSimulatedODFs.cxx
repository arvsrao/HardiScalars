/*
 *  Code for making an ODF image 
 *
 *
 *
 */
 
#include "itkSymRealSphericalHarmonicRep.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkVariableLengthVector.h"
#include "itkVectorContainer.h"
#include "itkImageFileWriter.h"
#include "itkNiftiImageIO.h"
//#include "SimulatedData.h"
#include "SimulatedCrosses.h"
#include "sphere1000points.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVariableLengthVectorCastImageFilter.h"
#include "string.h"

unsigned const short Dimension = 2;
unsigned const short NOrder = 8; 
	
	//Pixel Types
	typedef itk::SymRealSphericalHarmonicRep<double, NOrder>              PixelType;
	typedef itk::Image<PixelType, Dimension>                              ImageType;
	
	//IO types

	typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
	typedef itk::VectorImage<ImageType::PixelType::ComponentType, ImageType::ImageDimension>
                                                                                                   VectorImageType;
        typedef itk::VariableLengthVectorCastImageFilter<ImageType,VectorImageType>                CasterType;

        
        typedef itk::ImageFileWriter< VectorImageType >                                            WriterType;
         
int main(int argc, char* argv[]){
		
const char* outputfilename = argv[1];
	
	
	/*
                 * Convert resultDWI to a vector image
                 */
           
		ImageType::IndexType start; 
		ImageType::RegionType region;
                start[0] = 0;
	       	start[1] = 0;
	        		
		ImageType::SizeType size; 
		size[0] = 10;
		size[1] = 10;
                //size[2] = 1;
                region.SetSize( size );
                region.SetIndex( start );
                              
                ImageType::Pointer image = ImageType::New();
                image->SetRegions( region );
		image->Allocate();
          		
                //Declare an image iterator
		IteratorType iterator( image, region);
                iterator.GoToBegin(); 	       

               //Create PixelType
               PixelType     temp;
               
while(!iterator.IsAtEnd())
 {   
    //Fill temp pixel with RSH coefficents from simulated data file
     for(int i=0; i <45 ; ++i)
     {
        temp[i] = shCoef[iterator.GetIndex()[0] + iterator.GetIndex()[1]*size[0]][i];
     }
    std::cout<<" is placed at image index " <<  iterator.GetIndex()[0] + iterator.GetIndex()[1]*size[0] <<"\n";

    image->SetPixel(iterator.GetIndex(), temp);
    ++iterator;    
 }


/* 
 *
 * Cast OdfImageType to a VectorImage, which ITK can write to a NIFTI file
 *
 */

    //instantiate writer
  WriterType::Pointer writer = WriterType::New();
  writer->SetImageIO( itk::NiftiImageIO::New() );
  writer->SetFileName( outputfilename  );

  //Casting is done here
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( image  );
 
  printf("Casting ODF image\n");
  caster->Update();
              
try
 {
    writer->SetInput( caster->GetOutput() );
    printf("Writing ODF image\n");
    writer->Update();
 }
catch( itk::ExceptionObject & err )
 {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
 }
        
}//end main
	

	
        
