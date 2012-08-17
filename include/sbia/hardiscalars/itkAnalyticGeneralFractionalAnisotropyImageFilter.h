#ifndef __itkAnalyticGeneralFractionalAnisotropyImageFilter_h
#define __itkAnalyticGeneralFractionalAnisotropyImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix.h"
#include "math.h"

namespace itk
{
/** \class GeneralFractionalAnisotropyImageFilter
* \brief Implements an operator for normalizing images of real spherical
* harmonic pixels.
*
* \warning No numeric overflow checking is performed in this filter.
* \warning Pixels must implement a normalize method
*/

namespace Functor{  

template< class TInput, class TOutput>
class AnalyticGeneralFractionalAnisotropy{
public:

typedef vnl_vector_fixed<double, TInput::Dimension> OdfRepVectorType;

AnalyticGeneralFractionalAnisotropy() {}
virtual ~AnalyticGeneralFractionalAnisotropy() {}

inline TOutput operator() ( const TInput & A) 
{

OdfRepVectorType vec;
double lOne=0; 

//cast A to a vnl_vector_fixed
for(int i=0; i < TInput::Dimension ; ++i)
{
    vec.put(i , A[i]); 
}

if(vec.one_norm() ==0)
   return static_cast<TOutput>( 0 );
else
{
   //vec = vec / vec.two_norm();	
   lOne = LOneNormOfODF(vec);
   vec = vec/ lOne;
   //std::cout<< lOne  << std::endl;
    return static_cast<TOutput>( vcl_sqrt( dot_product(vec, vec) - .25 / vnl_math::pi ) );
  //return static_cast<TOutput>( vcl_sqrt( 1 - ( lOne *lOne /( 4 * vnl_math::pi * dot_product(vec,vec) ) ) ) );
  //return static_cast<TOutput>(lOne);
}
}

bool operator== (const AnalyticGeneralFractionalAnisotropy&) const
{
return true;
}

bool operator!= (const AnalyticGeneralFractionalAnisotropy&) const
{
return false;
}

private:

double factorial(double start, double end) //, int increment)
{

double prod=1;
if(end <=1)
{
   return 1;
}
else
{ 
   for(double i=start ; i <= end ; ++i)
   {
   	prod *= i;
   }
   return prod;
}
}//end factorial

double binomial(double a, double b)
{
    return factorial(1,a) / ( factorial(1,b) * factorial(1 , a-b));
}//end binomial 

double LOneNormOfODF(OdfRepVectorType vec)
{

double sum = 0;
double l_term = 0;
int index[5] = {0,3,10,21,36};

for(int l=0; l <= TInput::MaxOrder; l+=2)
{
      //if (l % 2 == 0) 
      for(int k=0; k <= l/2 ; ++k)
      {
        l_term += vcl_pow(-1, (double)k) * binomial(l, k) * factorial(l-2*k+2, 2*l-2*k);
      }//k-loop
      
      std::cout<<  l_term << std::endl;
      sum += ( l_term * vec.get(index[l/2]) ) / ( vcl_pow(2 ,(double)l) * factorial(1,l) ); 
      l_term = 0;
}//l-loop

return sum  * (double)2 * vnl_math::pi;
}//end function

};//end class
}//end functor

template <class TInputImage, class TOutputImage>
class ITK_EXPORT AnalyticGeneralFractionalAnisotropyImageFilter:
public

UnaryFunctorImageFilter<TInputImage,TOutputImage, 
Functor::AnalyticGeneralFractionalAnisotropy<typename TInputImage::PixelType,  typename TOutputImage::PixelType> >{
public:

/** Standard class typedefs. */

typedef AnalyticGeneralFractionalAnisotropyImageFilter  Self;
typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, Functor::AnalyticGeneralFractionalAnisotropy<typename TInputImage::PixelType,  
typename TInputImage::PixelType> >          Superclass;

typedef SmartPointer<Self>   Pointer;
typedef SmartPointer<const Self>  ConstPointer;

/** Method for creation through the object factory. */

itkNewMacro(Self);

/** Runtime information support. */

itkTypeMacro(AnalyticGeneralFractionalAnisotropyImageFilter, UnaryFunctorImageFilter);

protected:

AnalyticGeneralFractionalAnisotropyImageFilter() {}
virtual ~AnalyticGeneralFractionalAnisotropyImageFilter() {}

private:

AnalyticGeneralFractionalAnisotropyImageFilter(const Self&); //purposely not implemented
void operator=(const Self&); //purposely not implemented

};
} // end namespace itk

#endif
