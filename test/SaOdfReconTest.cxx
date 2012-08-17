#include <iostream>
#include "stdio.h"
#include "itkSymRealSphericalHarmonicRep.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkFixedArray.h"
#include "itkCastImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImage.h"
#include "itkSolidAngleOdfReconImageFilter.h"
#include "signal_mat.h"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_math.h"


/*double m_Delta = .001; 

double SmoothingFunction( double x )
{
    
    if( x < 0)
        return m_Delta / double(2);
    else if( x < m_Delta & x >=0)
        return m_Delta/double(2) + x*x/double(2*m_Delta);
    else if( x >=m_Delta & x < 1 - m_Delta )
        return x; 
    else if( x >= (1 - m_Delta) & x <1 )
        return 1 - m_Delta/double(2) - (1-x)*(1-x)/double(2*m_Delta) ;
    else 
        return 1 - m_Delta/double(2);   
}
*/
 
template <typename FilterType>
void addGradImage(FilterType *filter)
{
  const unsigned int numberOfGradientImages = 65; // The bbl set!!!

  typedef typename FilterType::GradientImagesType          GradImageType;
  typedef typename GradImageType::Pointer                  GradientImagePointer;
  typedef typename GradImageType::RegionType               GradientRegionType;
  typedef typename GradientRegionType::IndexType           GradientIndexType;
  typedef typename GradientRegionType::SizeType            GradientSizeType;
  typedef typename GradImageType::PixelType                GradPixelType;

  GradientImagePointer gradientImage = GradImageType::New();
  gradientImage->SetVectorLength(numberOfGradientImages);
  
  GradientSizeType  sizeGradientImage  = {{ 1, 1, 3}};
  GradientIndexType indexGradientImage = {{ 0, 0, 0}};
  GradientRegionType     regionGradientImage;
  regionGradientImage.SetSize(  sizeGradientImage );
  regionGradientImage.SetIndex( indexGradientImage);
  gradientImage->SetRegions( regionGradientImage );
  gradientImage->Allocate();

  itk::ImageRegionIteratorWithIndex< GradImageType > git(
        gradientImage, regionGradientImage );

  git.GoToBegin();

  //Create gradient pixel
  GradPixelType pix(numberOfGradientImages);
  
  // lets fill in the signal where we want to...

  //Odf1 in {0,0,0}
  indexGradientImage[0] = 0;
  indexGradientImage[1] = 0;
  indexGradientImage[2] = 0;
  pix.Fill(0);
  for(unsigned int i=0;i<numberOfGradientImages;++i)
  {
      pix[i] = signals[i][0];
  }
  gradientImage->SetPixel(indexGradientImage,pix);
    
    //Odf0 in {0,0,1}
    indexGradientImage[0] = 0;
    indexGradientImage[1] = 0;
    indexGradientImage[2] = 1;  
    pix.Fill(0);
    for(unsigned int i=0;i<numberOfGradientImages;++i)
    {
        pix[i] = signals[i][1] / signals[0][1];
        //  std::cout<< vcl_log( -1 * vcl_log( SmoothingFunction( pix[i] ))) << "\n";   
    }
    gradientImage->SetPixel(indexGradientImage,pix);
    
    //std::cout<< gradientImage->GetPixel(indexGradientImage);
    
    /*  for(unsigned int i=0; i< numberOfGradientImages; ++i)
     {
     std::cout<< signals[i][0]<< "\n";
     }
     */  

  ///Odf2 in {0,0,2}
  indexGradientImage[0] = 0;
  indexGradientImage[1] = 0;
  indexGradientImage[2] = 2;
  pix.Fill(0);
    for (unsigned int i=0;i<numberOfGradientImages;++i)
    {
       pix[i] = signals[i][2];
    }
    gradientImage->SetPixel(indexGradientImage,pix);
    
  //Set up Gradient Contatiner...
  typedef typename FilterType::GradientDirectionContainerType
                                          GradientDirectionContainerType;

  typedef typename FilterType::GradientDirectionType
                                          GradientDirectionType;
  
  typename GradientDirectionContainerType::Pointer gradContainer = GradientDirectionContainerType::New();

  GradientDirectionType dir;
  double  gradientDirections[numberOfGradientImages][3] =
  {
    {0.000000,0.000000,0.000000},
    {1.000000,0.000000,0.000000},
    {0.000000,1.000000,0.000000},
    {-0.026007,0.649170,0.760199},
    {0.591136,-0.766176,0.252058},
    {-0.236071,-0.524158,0.818247},
    {-0.893021,-0.259006,0.368008},
    {0.796184,0.129030,0.591137},
    {0.233964,0.929855,0.283956},
    {0.935686,0.139953,0.323891},
    {0.505827,-0.844710,-0.174940},
    {0.346220,-0.847539,-0.402256},
    {0.456968,-0.630956,-0.626956},
    {-0.486997,-0.388997,0.781995},
    {-0.617845,0.672831,0.406898},
    {-0.576984,-0.104997,-0.809978},
    {-0.826695,-0.520808,0.212921},
    {0.893712,-0.039987,-0.446856},
    {0.290101,-0.541189,-0.789276},
    {0.115951,-0.962591,-0.244896},
    {-0.800182,0.403092,-0.444101},
    {0.513981,0.839970,0.173994},
    {-0.788548,0.152912,-0.595659},
    {0.949280,-0.233069,0.211062},
    {0.232964,0.782880,0.576911},
    {-0.020999,-0.187990,-0.981946},
    {0.216932,-0.955701,0.198938},
    {0.774003,-0.604002,0.190001},
    {-0.160928,0.355840,0.920587},
    {-0.147035,0.731173,-0.666158},
    {0.888141,0.417066,0.193031},
    {-0.561971,0.231988,-0.793959},
    {-0.380809,0.142928,0.913541},
    {-0.306000,-0.199000,-0.931001},
    {-0.332086,-0.130034,0.934243},
    {-0.963226,-0.265062,0.044010},
    {-0.959501,0.205107,0.193101},
    {0.452965,-0.888932,0.067995},
    {-0.773133,0.628108,0.088015},
    {0.709082,0.408047,0.575066},
    {-0.692769,0.023992,0.720760},
    {0.681659,0.528735,-0.505747},
    {-0.141995,-0.724976,0.673978},
    {-0.740168,0.388088,0.549125},
    {-0.103006,0.822044,0.560030},
    {0.584037,-0.596038,0.551035},
    {-0.088008,-0.335031,0.938088},
    {-0.552263,-0.792377,0.259123},
    {0.838158,-0.458086,-0.296056},
    {0.362995,-0.560993,0.743990},
    {-0.184062,0.392133,-0.901306},
    {-0.720938,-0.692941,0.008999},
    {0.433101,0.682159,-0.589137},
    {0.502114,0.690157,0.521119},
    {-0.170944,-0.508833,-0.843722},
    {0.462968,0.422971,0.778946},
    {0.385030,-0.809064,0.444035},
    {-0.713102,-0.247035,0.656094},
    {0.259923,0.884737,-0.386885},
    {0.001000,0.077002,-0.997030},
    {0.037002,-0.902057,0.430027},
    {0.570320,-0.303170,-0.763428},
    {-0.282105,0.145054,-0.948354},
    {0.721098,0.608082,0.332045},
    {0.266985,0.959945,-0.084995}
  };
  
  for (unsigned int g = 0; g<numberOfGradientImages;++g)
  {
    dir[0] = gradientDirections[g][0];
    dir[1] = gradientDirections[g][1];
    dir[2] = gradientDirections[g][2];
    gradContainer->InsertElement(g,dir);
  }

  filter->SetGradientImage( gradContainer, gradientImage );
 /*
    //Print out signals
    //GradientDirectionContainerType  gradDirections  = reconFilter->GetGradientDirectionContainer(); 

    for(unsigned int i = 0; i < 3; ++i)
    {
        indexGradientImage[0] = 0;
        indexGradientImage[1] = 0;
        indexGradientImage[2] = i;
        std::cout<< gradientImage->GetPixel(indexGradientImage) << "\n\n";
    } 
  */  
}

template <typename FilterType>
void addMask(FilterType *filter)
{
  //We need to make a mask.
  typedef itk::Image< unsigned char , 3 >   ImageMaskType;
  ImageMaskType::Pointer maskImage = ImageMaskType::New();
  typedef ImageMaskType::RegionType MaskRegionType;
  typedef MaskRegionType::IndexType MaskIndexType;
  typedef MaskRegionType::SizeType  MaskSizeType;
  MaskSizeType  sizeMaskImage  = {{ 1, 1, 3 }};
  MaskIndexType indexMaskImage = {{ 0, 0, 0 }};
  MaskRegionType     regionMaskImage;
  regionMaskImage.SetSize(  sizeMaskImage );
  regionMaskImage.SetIndex( indexMaskImage);
  maskImage->SetRegions( regionMaskImage );
  maskImage->Allocate();
  maskImage->FillBuffer( static_cast< ImageMaskType::InternalPixelType >( 1 ) );

  typedef itk::ImageMaskSpatialObject< 3 >   MaskType;
  MaskType::Pointer  spatialObjectMask = MaskType::New();

  spatialObjectMask->SetImage( maskImage );
  filter->SetImageMask( spatialObjectMask );
    
  /*  //Testing mask creation
    MaskIndexType index;
    for(unsigned int j = 0; j<3 ;++j) 
    {
        index[0]=0; index[1]=0; index[2]=j;
        unsigned char a = maskImage->GetPixel(index);
        std::cout<< "value at index " << index << " is: "<< int(a) << "\n";
    }
    */
}

int SolidAngleOdfReconImTest1(int,char*[])
{
  //Just try to instantiate the filter...

  typedef unsigned int    ReferencePixelType;
  typedef unsigned int    GradientPixelType;
  typedef double          FopdfPrecisionType;
  typedef itk::SymRealSphericalHarmonicRep< FopdfPrecisionType, 4 > FopdfPixelType;

  typedef itk::SolidAngleOdfReconImageFilter<
    GradientPixelType,
    FopdfPixelType , 3 >      FopdfReconImageFilterType;

  try
  {
    FopdfReconImageFilterType::Pointer symReconFilter =
    FopdfReconImageFilterType::New();
      
    //Testing mask creation 
    std::cout<< "add the mask here" << "\n";
    addMask<FopdfReconImageFilterType>( symReconFilter );  

    //std::cerr << symReconFilter << std::endl;
    std::cout << "Success in instantiating filter" <<"\n";
    return EXIT_SUCCESS;
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

}//end first test

int SolidAngleOdfReconImTest2(int, char*[])
{

  //Test with laplace Beltrami smoothing
  
    typedef float           PrecisionType;

    /*Coefficients determined by Solid Angle Reconstruction (Sapiro, et al.) with Laplacian smoothing \delta = .01. 
      In particular:
             
                     shCoef = P*L*(B^{T}*B + .01*L^2)*BT * signal;
    */
    
     //Single diffusion tensor (1,1,10)  [ (x,y,z)]
    PrecisionType c1[45] = {0.2820947918,5.22875969797987e-05,0.000115206489072933,0.140960103745836,6.44024778232704e-05,2.90429209230895e-05,-3.90019909128842e-05,0.00063636818651428,0.000572455813324666,-0.000145273651185882,0.0302297067211656,-0.000126512740219173,-0.000128954410809254,-2.9127795696424e-05,0.00109322054617558,-0.000318046177058288,9.08642047262129e-05,0.000343898224919089,-0.000130399304871703,-5.57742958608139e-05,-0.000119991458650256,0.00383181705871875,0.000150988410992905,-0.000177669985265487,0.000152012818756404,0.00086399442448334,-0.000267531445562241,-0.000515945784246076,0.000686835952905499,-4.69611547206463e-05,-0.000288653062137111,0.000288384550273585,0.000150454596112847,0.000194032800281664,8.61748025294835e-05,0.000158366505873211,0.000211338838455017,0.000245450430028945,-0.000106393807320407,6.88933358347532e-05,0.000275027726476034,-0.000376891884660752,-0.000344178112160947,9.53017906860919e-05,-0.000230765782775903};

    //Single diffusion tensor (10,1,1)  [ (x,y,z)]
    PrecisionType c2[45] = {0.2820947918,0.122077516016263,-1.915342048906e-05,-0.071120715177173,6.1023733927947e-05,-0.000116509009758317,0.0228395751280926,-0.000214393044482409,-0.0174643664961144,-6.94410773494753e-06,0.0119437348611448,0.000310224781971073,-0.000252486515758973,-0.000272306542966044,-0.000410853202938509,0.00275691235539262,-0.000257016512104131,-0.00193137332474592,2.45437882368362e-06,0.0019543899388685,7.76889288676007e-05,-0.00131505030848115,0.000144660206539456,9.08898231163548e-05,0.00052438550271467,-7.52463445443056e-05,-2.02893336599579e-05,0.000320920114487901,0.000215519457304092,9.6451646392254e-05,-0.000140307582868918,0.000334002646533914,0.000132465722139676,-0.000203547792609675,-8.9033966009839e-05,-0.000119938239596427,0.000422845990234785,-0.000425298849878213,5.59060126716811e-05,-0.000422850210050815,-0.000184325874106101,0.000188161298810032,-9.47092192376993e-05,6.17605328509237e-05,0.00019963118608823};

    //Two diffusion tensors; sum of (1,1,10) & (10,1,1) with 90 degree difference between lobes.
    PrecisionType c3[45] = {0.2820947918,0.0564170002444126,0.00011450643378101,0.0310229845418324,5.88146107416707e-05,-0.000477243594106525,0.0257030062037759,0.000681300865192394,-0.0419237153527801,-0.000494461141561191,0.0684418087306622,0.000201545185532586,-0.000977288571441052,-0.000544921431225896,0.00078927453912132,0.00391185517026058,-0.000321538268210365,-0.00365419958165407,0.000134791223903194,0.00316139059269021,-0.000474911145655481,0.00449265896410093,0.000426836435631937,-0.000291010904805993,0.0010825047523546,0.00139878313524075,-0.000611446664734754,-0.000581019447971569,0.00136546627616936,7.98437896282088e-05,-0.000834164167614779,0.000937121154892438,-0.000370879197598778,-0.000143520508908867,0.00117397657216018,0.000119178521762813,-0.00017478852992334,-0.000490080936673768,-7.63237227915524e-05,-0.000415742788224926,0.000366209844839279,-0.000185062241577529,-0.000650600798307806,0.000338800486573612,0.0002441294680779};

  typedef float    GradientPixelType;
  typedef itk::SymRealSphericalHarmonicRep< PrecisionType, 8 >
                                                          PixelType;

  typedef itk::SolidAngleOdfReconImageFilter< GradientPixelType, PixelType,  3 >   
                                                                    ReconImageFilterType;
    
  typedef ReconImageFilterType::GradientDirectionContainerType     GradientDirectionContainerType;  


  ReconImageFilterType::Pointer reconFilter = ReconImageFilterType::New();
  reconFilter->SetBeltramiLambda(.01);
    reconFilter->SetDelta(.001);
    
  //Basic parameters...
  const double precision = 0.015;

  //Attach the grads and signal!
  addGradImage<ReconImageFilterType>(reconFilter);
  
  //Attach a mask 
   std::cout<< "add the mask here" << "\n";
   addMask<ReconImageFilterType>(reconFilter);
  
  //Make 3 RSHs to check....
  typedef ReconImageFilterType::OutputPixelType OutputPixelType;

  OutputPixelType pix1(c1);
  OutputPixelType pix2(c2);
  OutputPixelType pix3(c3);

  // Also see if vnl_svd is thread safe now...
  //std::cout << std::endl << "This filter is using " <<
  //reconFilter->GetNumberOfThreads() << " threads " << std::endl;
    
  try
  {
    reconFilter->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught on filter update!" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }

  typedef ReconImageFilterType::OutputImageType           OutputImageType;
  typedef OutputImageType::IndexType                      OutputImageIndexType;
  
  OutputImageType::Pointer outputImage = reconFilter->GetOutput();
 ///OK lets check the pixels we set explicitly..
  OutputImageIndexType    outputImageIndex;

  bool passed = true;

  OutputPixelType testPixel;
  
   /* std::cout<< "\n\nThese are the Coefficients"<<"\n";
    for (unsigned int i=0; i<3; ++i) 
     {
        outputImageIndex[0] = 0; outputImageIndex[1] = 0; outputImageIndex[2] = i; 
        testPixel = outputImage->GetPixel( outputImageIndex );
        std::cout<< testPixel  <<"\n";
     }  
   */
    
    float max_diff[3] = {0,0,0};   
    
  outputImageIndex[0] = 0; outputImageIndex[1] = 0; outputImageIndex[2] = 0;
  testPixel = outputImage->GetPixel(outputImageIndex);
  std::cout << "testing Pixel1" << std::endl;
  for (unsigned int i=0;i<OutputPixelType::Dimension;i++)
  {
    fprintf(stderr,"index %d : %f - %f = %f\n",i,testPixel[i],pix1[i],fabs(testPixel[i] - pix1[i]));
    
      if( max_diff[0] < fabs(testPixel[i] - pix2[i]) )
          max_diff[0] = fabs( testPixel[i] - pix2[i] );
  
      
    if (fabs(testPixel[i] - pix1[i]) > precision)
      passed = false;
  }

  outputImageIndex[0] = 0; outputImageIndex[1] = 0; outputImageIndex[2] = 1;
  testPixel = outputImage->GetPixel(outputImageIndex);
   
  std::cout << "testing Pixel2" << std::endl;
  for (unsigned int i=0;i<OutputPixelType::Dimension;i++)
  {
      fprintf(stderr,"index %d : %f - %f = %f\n",i,testPixel[i],pix2[i],fabs(testPixel[i] - pix2[i]));
      
      if( max_diff[1] < fabs(testPixel[i] - pix2[i]) )
          max_diff[1] = fabs( testPixel[i] - pix2[i] );
      
      if(fabs(testPixel[i] - pix2[i]) > precision)
      passed = false;
  }

  outputImageIndex[0] = 0; outputImageIndex[1] = 0; outputImageIndex[2] = 2;
  testPixel = outputImage->GetPixel(outputImageIndex);
  std::cout << "testing Pixel3" << std::endl;
  for (unsigned int i=0;i<OutputPixelType::Dimension;i++)
  {
    fprintf(stderr,"index %d : %f - %f = %f\n",i,testPixel[i],pix3[i],fabs(testPixel[i] - pix3[i]));
      
      if( max_diff[2] < fabs(testPixel[i] - pix2[i]) )
          max_diff[2] = fabs( testPixel[i] - pix2[i] );
  
    if (fabs(testPixel[i] - pix3[i]) > precision)
      passed = false;
  }


     FILE *fd;
    fd = fopen("/Users/arvind/Desktop/Test_Odf_Coefficents.h", "w");
    
    //Write to File Coefficients from testPixel
    for(unsigned int j =0; j<=2; ++j)
    {
        outputImageIndex[0] = 0; outputImageIndex[1] = 0; outputImageIndex[2] = j;
        testPixel = outputImage->GetPixel(outputImageIndex);
        
        fprintf(fd, "{ ");
        for(unsigned int i =0; i < OutputPixelType::Dimension ; ++i)
        { 
            if(i==OutputPixelType::Dimension-1)
                fprintf(fd, "%.15g},\n\n\n\n ", testPixel[i]);
            else 
                fprintf(fd, "%.15g, ", testPixel[i] ); 
        }
     /*   
        fprintf(fd, "[ ");
        for(unsigned int i =0; i < 65 ; ++i)
        { 
            if(i==65-1)
                fprintf(fd, "%.15g],\n\n\n\n ", signals[i][1]);
            else 
                fprintf(fd, "%.15g, ", signals[i][1] ); 
        }
      */
    }
    fclose(fd);  


    std::cout<< "\n\nlambda value is: "<< reconFilter->GetBeltramiLambda()<<"\n";
    std::cout<< "Precision is:  " << precision<<"\n";
    std::cout<< " Delta value is: "<< reconFilter->GetDelta()<<"\n";
    std::cout<< "\n\nThe maximum difference for pixel one is: "<< max_diff[0] <<"\n";
    std::cout<< "The maximum difference for pixel two is: "<< max_diff[1] <<"\n";
    std::cout<< "The maximum difference for pixel three is: "<< max_diff[2] <<"\n";
    
  if (!passed){
    std::cerr << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
  }

  std::cerr << "[PASSED]" << std::endl;
  return EXIT_SUCCESS;
}//end second test


   
int main( int argc, char * argv[])
{
  //if ( SolidAngleOdfReconImTest1(argc,argv) == EXIT_FAILURE )
  //  return EXIT_FAILURE;

  if ( SolidAngleOdfReconImTest2(argc,argv) == EXIT_FAILURE )
    return EXIT_FAILURE;

  return EXIT_SUCCESS;

}
