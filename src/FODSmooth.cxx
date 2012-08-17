/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: FODSmooth.cxx,v $
  Language:  C++
  Date:      $Date: 2012-06-22 18:20:38 -0400 (Fri, 22 Jun 2012) $
  Version:   $Revision: 33 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkAnisotropicDiffusionImageFilterwithPowersGradient.h"
#include "itkImage.h"
#include "itkVariableLengthVectorCastImageFilter.h"
#include "itkVectorCastImageFilter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkJoinImageFilter.h"
#include "itkImageToVectorImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkNiftiImageIO.h"

template <typename TPixelType, unsigned int TVectorDimension,
        unsigned int TPowersDimension >
int runSmoothingwithBothDims(const char * inputFilename, const char * powersFilename, const char * outputFilename, unsigned int smoothingIter, float conductanceParam )
{

  const   unsigned int        Dimension = 3;
  typedef itk::VectorImage< TPixelType, Dimension >    VectorImageType;
  typedef itk::Image<TPixelType,Dimension> ImageType;
    

  typedef itk::ImageFileReader< VectorImageType >  ReaderType;
  typedef itk::ImageFileWriter< VectorImageType >  WriterType;
  
  typename ReaderType::Pointer readerInput = ReaderType::New();
  typename ReaderType::Pointer readerPowers = ReaderType::New();
  typename WriterType::Pointer writer = WriterType::New();
  
   std::cout << "image dimension" << TVectorDimension << "powers dimension" << TPowersDimension << std::endl;

  readerInput->SetFileName( inputFilename  );
  readerPowers->SetFileName( powersFilename  );
    try 
    { 
    readerInput->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 
      try 
    { 
    readerPowers->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    }  
  
   
  typedef itk::Vector<TPixelType,TVectorDimension> InputVectorType;
  typedef itk::Vector<TPixelType,TPowersDimension> PowersVectorType;
  typedef itk::Image< InputVectorType, Dimension >    InputImageVectorType;
  typedef itk::Image< PowersVectorType, Dimension >    PowersImageVectorType;
  
  
  
  typedef itk::VectorCastImageFilter<VectorImageType,InputImageVectorType> VecFilterType2;
  typename VecFilterType2::Pointer vecFilter2 = VecFilterType2::New();
  vecFilter2->SetInput(readerInput->GetOutput());
 
 
  typedef itk::VectorCastImageFilter<VectorImageType,PowersImageVectorType> PowersVecFilterType2;
  typename PowersVecFilterType2::Pointer powersvecFilter2 = PowersVecFilterType2::New();
  powersvecFilter2->SetInput(readerPowers->GetOutput());
     

  typedef itk::JoinImageFilter<PowersImageVectorType,InputImageVectorType> JoinFilterType;
  typename JoinFilterType::Pointer joinFilter = JoinFilterType::New();
  joinFilter->SetInput1(powersvecFilter2->GetOutput());
  joinFilter->SetInput2(vecFilter2->GetOutput());
 
 
  std::cout<< "creating smoothing filter" << std::endl;
  
  typedef typename JoinFilterType::OutputImageType JoinFilterOutputImageType;
  
   typedef itk::AnisotropicDiffusionImageFilterwithPowersGradient<JoinFilterOutputImageType,JoinFilterOutputImageType> SmoothFilterType;
  typename SmoothFilterType::Pointer smoother = SmoothFilterType::New();
  smoother->SetDrivingComponents(TPowersDimension);
  smoother->SetNumberOfIterations(smoothingIter);
  smoother->SetConductanceParameter(conductanceParam);
  smoother->SetTimeStep(0.05);  //for .125 received a warning: That time should be < .0625
  smoother->SetInput(joinFilter->GetOutput());
  smoother->SetDebug(false);
 itk::SimpleFilterWatcher watcher(smoother,"smoother" ); 
 
   try 
    { 
    smoother->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 
  
    
  typedef itk::VariableLengthVectorCastImageFilter<JoinFilterOutputImageType,VectorImageType> VecFilterType;
  typename VecFilterType::Pointer vecFilter = VecFilterType::New();
  vecFilter->SetInput(smoother->GetOutput());
   
  typedef itk::VectorIndexSelectionCastImageFilter<VectorImageType,ImageType> IndexSelectionFilterType;
  typename IndexSelectionFilterType::Pointer indexFilter = IndexSelectionFilterType::New();
 indexFilter->SetInput(vecFilter->GetOutput() );
 
 typedef itk::ImageToVectorImageFilter<ImageType> ImVecFilterType;
  typename ImVecFilterType::Pointer imvecFilter = ImVecFilterType::New();
  typename ImageType::Pointer image = ImageType::New();
 
 for (unsigned int i= TPowersDimension; i<TVectorDimension+TPowersDimension ; i++)
 { 
 
   indexFilter->SetIndex(i);
   std::cout << "Reading Volume .... "<< i << std::endl;
   try 
    { 
    indexFilter->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 
		  
  image = indexFilter->GetOutput(); 
  imvecFilter->SetNthInput( i-TPowersDimension, image );   
  image->DisconnectPipeline();
  }

 try
   { 
  imvecFilter->Update();
   }
 catch (itk::ExceptionObject &ex)
   {
     std::cout << ex << std::endl;
   } 

    
   
  std::cout<< "writing it out" << std::endl;
  writer->SetImageIO( itk::NiftiImageIO::New() );
  writer->SetFileName( outputFilename );
  writer->SetInput( imvecFilter->GetOutput() );
  try 
    { 
    writer->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 
    
   
  return EXIT_SUCCESS;

}


template <typename TPixelType, unsigned int TVectorDimension >
int runSmoothingwithVectorDim(const char * inputFilename, const char * powersFilename, const char * outputFilename, unsigned int smoothingIter, float conductanceParam, unsigned int PowersDimension )
{
   switch (PowersDimension)
     {
//     case 1:
//       return runSmoothingwithBothDims<TPixelType,TVectorDimension, 1>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam);
//     	break;
     case 2:
       return runSmoothingwithBothDims<TPixelType,TVectorDimension, 2>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam);
     	break;
     case 3:
      return runSmoothingwithBothDims<TPixelType,TVectorDimension, 3>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam);
     	break;
     case 4:
       return runSmoothingwithBothDims<TPixelType,TVectorDimension, 4>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam);
     	break;
     case 5:
       return runSmoothingwithBothDims<TPixelType,TVectorDimension, 5>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam);
     	break;
     default:
       return runSmoothingwithBothDims<TPixelType,TVectorDimension, 5>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam);
     }  
    
     return EXIT_SUCCESS;
}


int main( int argc, char ** argv )
{
  // Verify the number of parameters in the command line
  if( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImageFile RshPowersFile outputImageFile " << std::endl;
    return EXIT_FAILURE;
    }

  const char * inputFilename  = argv[1];
  const char * powersFilename = argv[2];
  const char * outputFilename = argv[3];
  unsigned int smoothingIter=atoi(argv[4]);
  float conductanceParam=atof(argv[5]);  
  typedef float      PixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::VectorImage< PixelType, Dimension >    VectorImageType;
  typedef itk::ImageFileReader< VectorImageType >  ReaderType;
  ReaderType::Pointer readerInput = ReaderType::New();
    
  readerInput->SetFileName( inputFilename  );
   try 
    { 
    readerInput->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 
   
    
  ReaderType::Pointer readerPowers = ReaderType::New();
    
  readerPowers->SetFileName( powersFilename  );
   try 
    { 
    readerPowers->Update(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
    } 
    
  unsigned int PowersDimension=readerPowers->GetOutput()->GetVectorLength();    
  unsigned int VectorDimension=readerInput->GetOutput()->GetVectorLength();
     std::cout << "image dimension" << readerInput->GetOutput() << "powers dimension" << PowersDimension << std::endl;

     switch (VectorDimension)
     {
//     case 1:
//       return runSmoothingwithVectorDim<PixelType, 1>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam,PowersDimension);////
//     	break;
     case 2:
       return runSmoothingwithVectorDim<PixelType, 2>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam,PowersDimension);
     	break;
     case 3:
       return runSmoothingwithVectorDim<PixelType, 3>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam,PowersDimension);
     	break;
     case 4:
       return runSmoothingwithVectorDim<PixelType, 4>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam,PowersDimension);
     	break;
     case 5:
       return runSmoothingwithVectorDim<PixelType, 5>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam,PowersDimension);
     	break;
     case 6:
       return runSmoothingwithVectorDim<PixelType, 6>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam,PowersDimension);
     	break;
     case 15:
       return runSmoothingwithVectorDim<PixelType, 15>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam,PowersDimension);
     	break;
     case 28:
       return runSmoothingwithVectorDim<PixelType, 28>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam,PowersDimension);
     	break;
     case 45:
       return runSmoothingwithVectorDim<PixelType, 45>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam,PowersDimension);
     	break;
     default:
       return runSmoothingwithVectorDim<PixelType, 45>(inputFilename,powersFilename,outputFilename,smoothingIter,conductanceParam,PowersDimension);
     }
     
      
 
     return EXIT_SUCCESS;
}
