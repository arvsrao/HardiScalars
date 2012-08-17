#ifndef __itkCalculateOdfFromSamplesImageFilter_h
#define __itkCalculateOdfFromSamplesImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix.h"

namespace itk

{


/** \class CalculateOdfFromSamplesImageFilter
* \brief Implements an operator for normalizing images of real spherical
* harmonic pixels.
*
* 
* \warning No numeric overflow checking is performed in this filter.
* \warning Pixels must implement a normalize method
*/

namespace Functor {  

template< class TInput, class TOutput>
class CalculateOdfFromSamples{
	public:

typedef vnl_vector_fixed<double, TInput::Dimension> OdfRepVectorType;


CalculateOdfFromSamples(){

int N = TInput::PixelType::Size();
int n = vnl_sqrt(N);

typedef vnl_vector_fixed<double, const int N> VectorType;
typedef vnl_matrix<double, TOutput::Dimension, N> MatrixType; 

VectorType weight; 
MatrixType mat; 

ComputeOdfFromSamples(){

VectorType phi;
VectorType theta;
int M = TOutput::NOrder;

//create phi vector
for( int i = 0 ; i < n ){ 
   phi.put( i, i * 2 * M_PI / n );    //delta = \pi  / N 
}


//create theta vector
for(int i = 0; i < n ){
  theta.put( i, i * M_PI / n );     //delta = .5 * \pi /N
}


//create SH basis matrix
for( int i=0; i < TOutput::Dimension ; ++i){
   for( int j=0; j < n ; ++j){
      for( int k = 0; j < n ; ++k){
         mat.put(i, k +j*n, TOutput::Y(i, theta.get(j ), phi,get( k ) ));  
}}}

//create weight vector
for(int i=0; i< n; ++i){
  for(int j=1; j<= n; ++j){
     weight.put(j + n*i, ( phi.get(j) - phi.get(j-1) ) * ( vcl_cos(theta.get(i)) - vcl_cos(theta.get(i-1))) ); 
}}

}//end constructor



virtual ~CalculateOdfFromSamples(){}


inline TOutput operator() ( const TInput & A ) const
{

return static_cast<TOutput>(mat * ( ) );
}



bool operator== (const CalculateOdfFromSamples&) const
{

return true;
}

bool operator!= (const CalculateOdfFromSamples&) const
{

return false;
}


};  //end Class


template <class TInputImage, class TOutputImage, int N>
class ITK_EXPORT CalculateOdfFromSamplesImageFilter :
 public

UnaryFunctorImageFilter<TInputImage,TOutputImage, 

Functor::CalculateOdfFromSamples<typename TInputImage::PixelType,  typename TOutputImage::PixelType, N  > > 

{

public:

/** Standard class typedefs. */

typedef CalculateOdfFromSamplesImageFilter  Self;

typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 

Functor::CalculateOdfFromSamples<typename TInputImage::PixelType,  

typename TInputImage::PixelType, N  > >          Superclass;


typedef SmartPointer<Self>   Pointer;

typedef SmartPointer<const Self>  ConstPointer;


/** Method for creation through the object factory. */

itkNewMacro(Self);


/** Runtime information support. */

itkTypeMacro(CalculateOdfFromSamplesImageFilter, 

UnaryFunctorImageFilter);


protected:

CalculateOdfFromSamplesImageFilter() {}

virtual ~CalculateOdfFromSamplesImageFilter() {}


private:

CalculateOdfFromSamplesImageFilter(const Self&); //purposely not implemented

void operator=(const Self&); //purposely not implemented


};


} // end namespace itk

#endif
