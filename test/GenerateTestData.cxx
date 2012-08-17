/*
 *  GenerateTestDate.cxx
 *  
 *
 *  Created by Arvind  Rao on 11/29/10.
 *   
 *  The purpose of this file is generate an 2-D ODF image 
 *  to test HARDI scalar measures. 
 *
 */

#include "itkImage.h"
#include "itkVector.h"
#include "itkVectorImage.h"
#include "itkVariableLengthVector.h"
#include "itkVectorContainer.h"
#include "itkImageFileWriter.h"
#include "itkNiftiImageIO.h"
#include "itkSymRealSphericalHarmonicRep.h"
#include "vnl/vnl_vector_fixed.h"
#include "itkOdfReconImageFilter.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix_fixed.h"
#include "sphere1000points.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVariableLengthVectorCastImageFilter.h"
#include "string.h"

/****************************************************
	 Typedef's and Constants
****************************************************/
	
	unsigned const short Dimension = 3;
	unsigned const short NOrder = 8; 
	
	typedef vnl_vector_fixed<double, 1000> VectorType;	
        typedef vnl_matrix_fixed<double, 1000, 3> MatrixType; 
        //typedef itk::Vector<double, 1000> VecType; 

	//Pixel Types
	typedef itk::SymRealSphericalHarmonicRep<double, NOrder>              OdfPixelType;
	
	//Recon ODF filter
	typedef itk::OdfReconImageFilter<double, OdfPixelType>                OdfReconFilterType;
        typedef OdfReconFilterType::GradientImagesType                        DWI;
	typedef OdfReconFilterType::OutputImageType                           OdfImageType;

	//IO types
	//typedef itk::ImageFileWriter<DWI> TestType; 
	typedef itk::ImageRegionIteratorWithIndex< DWI > IteratorType;

	typedef itk::VectorImage<OdfImageType::PixelType::ComponentType, OdfImageType::ImageDimension>
                                                                                       VectorImageType;

        typedef itk::VariableLengthVectorCastImageFilter<OdfImageType,VectorImageType> CasterType;

   
        typedef itk::ImageFileWriter< VectorImageType >                                WriterType;
	

			
/****************************************************
                       MAIN
 ****************************************************/
	
	int main(int argc, char* argv[]){
		
	const char* outputfilename = argv[1];

		
	   /*
		TestType::Pointer writer = TestType::New();
		writer->SetImageIO( itk::NiftiImageIO::New() );
		writer->SetFileName( outputfilename );
   */
		
	        //std::cout << points1000X[1] << std::endl;

		/*
		 * Intially lets try to make a 4 by ten image 
		 * 
		 * Row #1: single lobe ODF's of different strengths and orientations. 
		 * Row #2: two lobe ODF's of different strengths and orientations 
			 * Row #3: two lobe ODF's of same strength and different orientations
		 * Row #4: three lobe ODF;s of different strengths and orientations
		 * 
		 */
		
		//cast points to a vnl_matrix
                MatrixType points; 
		VectorType theta;
		VectorType phi;

		for(int i =0; i <1000; ++i)
		{
                 points.put(i , 0 , points1000X[i]);
		 points.put(i , 1 , points1000Y[i]); 
		 points.put(i , 2 , points1000Z[i]); 
		}
	
         /*       //Fill theta and phi
		for(int i=0; i < 1000; ++i)
		{
		   theta.put( i, vcl_acos( points1000Z[i] ) );

		   if( (points1000Z[i]  != 1) & (points1000Z[i] != -1) )
		   {
		       if(( points1000X[i] ==0) & (points1000Y[i] < 0) )
			  phi.put(i, -1 * vnl_math::pi);
		       else if( (points1000X[i]==0) & (points1000Y[i] > 0) )
                          phi.put(i, vnl_math::pi);
                       else 
		          phi.put(i, vcl_atan( points1000Y[i] / points1000X[i]) ) ;
	           }
                   else 
		       phi.put(i, 0);
		}
*/
	        
		                	        
	        // Set of Peak/Mean vectors, row oriented.                 
		double directions_array[30] = {.99120, .05200, .12175, .72408, .27088, .63428, .18038, .38630, .90453, -.43222, .35416, .82929, -.87972, .18675, .43728, -.99120, -.05200, -.12175, -.72408, -.27088,-.63428,-.18038,-.38630,-.90453,.43222, -.35416, -.82929, .87972, -.18675, -.43728};
                /*
		 * (1,0,0) -- isotropic Odf
		 * (0,1,0) -- mean direction of (1, 0, 0) 
		 * (0,0,1) -- mean direction of (1, 1, 0 )
		 */   
               	double alpha_array[100] = { 
					    3,2,1,0,0,0,0,0,0,0,
                                            3,0,2,1,0,0,0,0,0,0,
                                            3,0,2,0,0,1,0,0,0,0,
                                            3,0,0,0,2,0,0,1,0,0,
                                            3,0,2,0,0,0,1,0,0,0,
					    3,0,0,2,0,0,0,0,1,0,
					    3,0,0,0,2,0,1,0,0,0,
					    3,0,0,0,0,0,0,2,0,1,
                                            3,0,0,1,0,0,0,0,0,2,

					    3,0,0,2,0,0,1,0,0,0,
                                        /*  0,2,0,1,0,0,0,0,0,0,
                                            0,2,0,0,1,0,0,0,0,0,
                                            0,2,0,0,0,1,0,0,0,0,
                                            0,2,0,0,0,0,1,0,0,0,
					    0,2,0,0,0,0,0,1,0,0,
					    0,2,0,0,0,0,0,0,1,0,
					    0,2,0,0,0,0,0,0,0,1,

					    0,0,2,1,0,0,0,0,0,0,
                                            0,0,2,0,1,0,0,0,0,0,
                                            0,0,2,0,0,1,0,0,0,0,
                                            0,0,2,0,0,0,1,0,0,0,
                                            0,0,2,0,0,0,0,1,0,0,
					    0,0,2,0,0,0,0,0,1,0,
					    0,0,2,0,0,0,0,0,0,1,

					    0,0,0,2,1,0,0,0,0,0,
                                            0,0,0,2,0,1,0,0,0,0,
                                            0,0,0,2,0,0,1,0,0,0,
                                            0,0,0,2,0,0,0,1,0,0,
                                            0,0,0,2,0,0,0,0,1,0,
					    0,0,0,2,0,0,0,0,0,1,

					    0,0,0,0,2,1,0,0,0,0,
                                            0,0,0,0,2,0,1,0,0,0,
                                            0,0,0,0,2,0,0,1,0,0,
                                            0,0,0,0,2,0,0,0,1,0,
                                            0,0,0,0,2,0,0,0,0,1,
					    
                                            0,0,0,0,0,2,1,0,0,0,
                                            0,0,0,0,0,2,0,1,0,0,
                                            0,0,0,0,0,2,0,0,1,0,
                                            0,0,0,0,0,2,0,0,0,1,
                                            
					    0,0,0,0,0,0,2,1,0,0,
                                            0,0,0,0,0,0,2,0,1,0,
                                            0,0,0,0,0,0,2,0,0,1,
                                            
					    0,0,0,0,0,0,0,2,1,0,
                                            0,0,0,0,0,0,0,2,0,1,

                                            0,0,0,0,0,0,0,0,2,1, */

            };
		
	        double test[100];
		vnl_matrix_fixed<double, 10, 10>               alpha;  
		vnl_matrix_fixed<double, 1000, 10>             Lobes;                
                vnl_matrix_fixed<double, 10, 3>                MatrixOfMeanDirections; 
		vnl_vector_fixed<double, 2>                    kappa(5, 5); //concentration parameter
	        
                //Perform transpose of alpha_array.		
                for(int i =0; i < 10; ++i)
		{
			for(int j=0; j<10 ; ++j)
			{
		         test[i*10 + j] = alpha_array[i + j*10]; 		
			}
		}

		alpha.copy_in( test );
		MatrixOfMeanDirections.copy_in( directions_array );		
		

               /*
                * Create lobe matrix 1000 x Rank(alpha). Then Lobes * alpha = desired DWI
                */
		
                double C = .3957123096;

		vnl_vector_fixed<double,3> scratch(1,1,1); 
		for(int j = 0; j < 10; ++j)
	         {
		    for(int i = 0; i <1000; ++i)
		     {
		       scratch.put(0, -1 * MatrixOfMeanDirections.get_row(j).get(0));
		       scratch.put(1, -1 * MatrixOfMeanDirections.get_row(j).get(1));
		       scratch.put(2, -1 * MatrixOfMeanDirections.get_row(j).get(2));

		       Lobes.put(i,j, vcl_exp( 10 *  vcl_cos( dot_product( points.get_row(i), MatrixOfMeanDirections.get_row(j))) * vcl_cos( dot_product( points.get_row(i),  MatrixOfMeanDirections.get_row(j))) )/C );
                      /* Lobes.put(i, j, vcl_exp( -1*kappa.get(1) * vcl_sin( dot_product( points.get_row(i), scratch))) / (4 * vnl_math::pi * vcl_sinh(kappa.get(1)))
        + vcl_exp(-1* kappa.get(1) * vcl_sin( dot_product( points.get_row(i), MatrixOfMeanDirections.get_row(j)))) / (4 * vnl_math::pi * vcl_sinh(kappa.get(1))) );
                       */
		     }
                 }
	       vnl_matrix_fixed<double, 1000, 10> resultDWI = Lobes * alpha;
               
		/*
                 * Convert resultDWI to a vector image
                 */
           
		DWI::IndexType start; 
		DWI::RegionType region;
                start[0] = 0;
	       	start[1] = 0;
	        start[2] = 0;	
		
		DWI::SizeType size; 
		size[0] = 10;
		size[1] = 1;
                size[2] = 1;
                region.SetSize( size );
                region.SetIndex( start );
                              
                DWI::Pointer image = DWI::New();
                image->SetRegions( region );
		image->SetVectorLength(1001);
		image->Allocate();
          		
                itk::VariableLengthVector<double> temp;
	        temp.SetSize(1001);

	        //Declare an image iterator
		IteratorType iterator( image, region);
                iterator.GoToBegin(); 	       

while(!iterator.IsAtEnd())
 {   

   
    temp[0] = 1; //Baseline B_0 image
    std::cout<< "ODF with parameters " << alpha.get_column(iterator.GetIndex()[0]) << " is placed at image index " <<  iterator.GetIndex() <<"\n";

        for(int k = 1; k <=1000; ++k)
     {
	    	temp[k] = resultDWI.get(k-1, iterator.GetIndex()[0] );
		     }
    image->SetPixel(iterator.GetIndex(), temp);
    ++iterator;    
 }

                //Run ODF Recon Filter 
		typedef itk::VectorContainer<unsigned int, vnl_vector_fixed<double,3> > GradientDirectionContainerType;
                GradientDirectionContainerType::Pointer container = GradientDirectionContainerType::New();
		GradientDirectionContainerType::Iterator it = container->Begin();
            
	        vnl_vector_fixed<double,3> zero_vector(0,0,0); //needed to indicate at least one baseline image 
                for(int i =0 ; i <= 1000 ; ++i, ++it)
		 {
	          if(i==0)
		   {
	            container->InsertElement(it.Index(), zero_vector);
		   }
                  else
		   {		  
	            container->InsertElement(it.Index(), points.get_row(i));
		   }

		 }

		OdfReconFilterType::Pointer filter = OdfReconFilterType::New();
		filter->SetGradientImage(container, image); 
	        printf("computing ODF image\n");
        	filter->Update();
  
//Iterate through image an calculate GV 


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
  caster->SetInput( filter->GetOutput() );
 
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

	
	/*
		  vec.SetElement(0, points.get_row(i).get(0));
		  vec.SetElement(1, points.get_row(i).get(1));
		  vec.SetElement(2, points.get_row(i).get(2));
	
                for(int i=0; i <3; ++i)
		{ 
	            index[0]=i;
                    for(int j=0; j <10; ++j)
		    {
	                index[1]=j;
		        std::cout<< index <<std::endl;

		        for(int k = 0; k < 1000; ++k)
		        {
	                   temp[k] = resultDWI.get_column(i*10+j).get(k);
                        }

	                image->SetPixel(index, temp);
	                std::cout<< temp <<std::endl;	
		    }
		 }
           
 // Set of Peak/Mean vectors, row oriented.                 
		double directions_array[30] = { -0.376012006313, 0.3660250812, 0.851258251673  , 
			                          0.944455444004, -0.0231777516681, -0.327821149591,
						  0.266510503623, -0.878411539144, -0.396692978709,
						  0.402824982446, -0.633766117991, 0.660357888727,
						  -0.667212336099, -0.60320143784, 0.437007693233,
						  0.672415969315, 0.405588245249, 0.619156635696,
						  -0.672393437203, -0.4055725778, -0.619191367626,
						  0.165587926348, 0.976645609674, -0.1369086986,
						  0.0987799684209, 0.260343406975, -0.960449805187,
						  -0.834957885122, 0.535526475997, -0.1267151276

					}


					double alpha_array[300] =    {1,0,0,0,0,0,0,0,0,0,
		                              0,1,0,0,0,0,0,0,0,0,
					      0,0,1,0,0,0,0,0,0,0,
					      0,0,0,1,0,0,0,0,0,0,
					      0,0,0,0,1,0,0,0,0,0,
					      0,0,0,0,0,1,0,0,0,0,
					      0,0,0,0,0,0,1,0,0,0,
					      0,0,0,0,0,0,0,1,0,0,
					      0,0,0,0,0,0,0,0,1,0,
					      0,0,0,0,0,0,0,0,0,1,

					      .5,.5,0,0,0,0,0,0,0,0,
					      .5,0,.5,0,0,0,0,0,0,0,
                                              .5,0,0,.5,0,0,0,0,0,0,
					      .5,0,0,0,.5,0,0,0,0,0,
					      .5,0,0,0,0,.5,0,0,0,0,
                                              .5,0,0,0,0,0,.5,0,0,0,
                                              .5,0,0,0,0,0,0,.5,0,0,
                                              .5,0,0,0,0,0,0,0,.5,0,
                                              .5,0,0,0,0,0,0,0,0,.5,
				
	                                      .1,.9,0,0,0,0,0,0,0,0,
					      .2,.8,0,0,0,0,0,0,0,0,
                                              .3,.7,0,0,0,0,0,0,0,0,
					      .4,.6,0,0,0,0,0,0,0,0,
					       0,.2,0,0,.6,.2,0,0,0,0,
                                              .5,0,0,.1,0,0,0,0,0,0,.4,
                                               0,0,.4,0,0,.1,0,.5,0,0,
                                               0,.1,0,.3,0,.2,0,.4,0,0,
                                               0,0,.3,0,.3,0,.4,0,0,0 };
						  
*/
