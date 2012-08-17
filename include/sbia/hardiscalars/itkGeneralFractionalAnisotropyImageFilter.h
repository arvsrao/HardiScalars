#ifndef __itkGeneralFractionalAnisotropyImageFilter_h
#define __itkGeneralFractionalAnisotropyImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkNumericTraits.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_matrix.h"
#include "iostream.h"
#include "sphere1000points.h"
#include "itkMacro.h"

namespace itk
{
  
/** \class GeneralFractionalAnisotropyImageFilter
 * \brief Implements an operator for GFA calculation
 * 
 *
 * 
 * \warning No numeric overflow checking is performed in this filter.
 * \warning Pixels must implement a normalize method
 */

namespace Functor {  
  
template< class TInput, class TOutput >
class GeneralFractionalAnisotropy{

public:

typedef TInput                                                   PixelType; 
typedef typename PixelType::ComponentType                        ComponentType;
typedef TOutput                                                  OutputPixelType;

/** Dimension of the RSH pixel and Number of sphere samples */
itkStaticConstMacro(RshDimension, unsigned int, PixelType::Dimension );
itkStaticConstMacro(nmOfPoints, unsigned int, 1000 );

typedef vnl_matrix< ComponentType >                                    MatrixType;
typedef vnl_vector_fixed< ComponentType, nmOfPoints>                   VectorType;
ComponentType                                                          epsilon;

GeneralFractionalAnisotropy()
{
  std::cout << "constructing the Filter" <<"\n";
  epsilon = .0001;
  mat.set_size(nmOfPoints, RshDimension);
 
  m_Norm = false;  //Default behavior          
  
  //Fill theta and phi vectors
    for(int i=0; i < nmOfPoints; ++i)
    {
         theta.put( i, vcl_acos( points1000Z[i] ) );

         if( (points1000Z[i] == 1) || (points1000Z[i] == -1) ) //if you're at a pole
         {
            phi.put(i, 0);
         }
         else 
         {
            phi.put(i, vcl_atan2( points1000Y[i], points1000X[i] ) );
         }      
    }


     //create mat matrix of sampled SH basis functions
     for(int i=0; i < nmOfPoints; ++i)
      {
         for(int j = 1; j <= RshDimension; ++j)
          { 
             mat.put(i, j-1, TInput::Y(j, theta.get(i), phi.get(i) ) );
          }     
      }
}//end constructor

virtual ~GeneralFractionalAnisotropy() {} //destructor
/*
 *The 'const' suffix key word of the operator below forbids the writing of memeber variables.  
 *So, if VectorType samples is declared outside of operator(), then samples can't be change.
 *
 */
inline TOutput operator()( const TInput& A ) 
{  

vnl_vector_fixed< ComponentType, RshDimension> input;
PixelType B=A;

   // 'cast' A to a vnl_vector_fixed
   // ...also normalize RSH vector 
  if(m_Norm) 
   { 
     B.Normalize();
     for(int i=0; i < RshDimension; ++i) input.put( i , B[i] );
   }
  else 
     for(int i=0; i < RshDimension; ++i) input.put( i , A[i] ); 
     
      
   if( input.one_norm() < epsilon )
       return static_cast<TOutput>( 0 );
   else
   {
       
   VectorType                          samples = mat * input; 
   ComponentType                               sstd, rms, mean;
     
   mean = samples.mean();
   rms  = dot_product(samples, samples);
   sstd = rms - float(2)*mean*samples.sum() + ( mean * mean * float(nmOfPoints) ); 

   sstd /= float(nmOfPoints -1);
   rms /= float(nmOfPoints);
  
#ifdef TESTING 
   std::cout<< samples.get(0) << "\n";
   std::cout<< "root mean square: "<< rms << "\n";
   std::cout<< "sample std : "<< sstd << "\n";
#endif
        
   return static_cast<TOutput>( vcl_sqrt( sstd / rms ) );
   }
}

bool operator== (const GeneralFractionalAnisotropy&) const
  {
    return true;
  }
bool operator!= (const GeneralFractionalAnisotropy&) const
  {
    return false;
  }

void NormalizationOn(){m_Norm = true;}
void NormalizationOff(){m_Norm = false;}
void GetNormalizationUse(){return m_Norm;}

private:

MatrixType mat; 
VectorType theta; 
VectorType phi;

bool m_Norm; //should the RSH functions be normalized. 

}; //end Class
} //end namespace Functor


template <class TInputImage, class TOutputImage>
class ITK_EXPORT GeneralFractionalAnisotropyImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                       Functor::GeneralFractionalAnisotropy< typename TInputImage::PixelType,  typename TOutputImage::PixelType  > > 
{
public:
  /** Standard class typedefs. */
  typedef GeneralFractionalAnisotropyImageFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 
              Functor::GeneralFractionalAnisotropy<typename TInputImage::PixelType,  
              typename TInputImage::PixelType  > >          Superclass;
              
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set Normalization Flag; Defines On/Off Methods **/
  void NormalizationOn(){this->GetFunctor().NormalizationOn();}
  void NormalizationOff(){this->GetFunctor().NormalizationOff();}
  void GetNormalizationUse(){this->GetFunctor().GetNormalizationUse();} 


  /** Runtime information support. */
  itkTypeMacro(GeneralFractionalAnisotropyImageFilter, 
               UnaryFunctorImageFilter);

protected:
  GeneralFractionalAnisotropyImageFilter(){}
  virtual ~GeneralFractionalAnisotropyImageFilter() {}

private:
  GeneralFractionalAnisotropyImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};
} // end namespace itk


#endif
