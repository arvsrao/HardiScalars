
#ifndef __itkSolidAngleOdfReconImageFilter_txx
#define __itkSolidAngleOdfReconImageFilter_txx

#include "itkSolidAngleOdfReconImageFilter.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_math.h"
#include "itkReplaceSpecialFunctions.h"

namespace itk 
{

template< 
          class TGradientImagePixelType,
          class TOutputPixelType,
          unsigned int TImageDimension,
          class TResidualPercisionType
        >
SolidAngleOdfReconImageFilter< TGradientImagePixelType, TOutputPixelType, TImageDimension, TResidualPercisionType >
::SolidAngleOdfReconImageFilter()
{
  this->AllowSuperResolutionOff();
  m_Delta = .001;
}

template< 
          class TGradientImagePixelType,
          class TOutputPixelType,
          unsigned int TImageDimension,
          class TResidualPercisionType
        >
vnl_vector< double >
SolidAngleOdfReconImageFilter<  TGradientImagePixelType, TOutputPixelType, TImageDimension, TResidualPercisionType >
::ComputeCoeffsFromSignal( vnl_vector< double > signal, ResidualPixelType& residual )
{
  //call the superclass method to compute the RSH coeffs of the signal!
  
  vnl_vector< double > log_signal( signal.size() );
  //log_signal.set_size( signal.size());
  
  for(unsigned int i = 0; i < signal.size(); ++i)
   {
     log_signal.put(i, vcl_log( -1 * vcl_log( SmoothingFunction(signal.get(i)))) );  
   }  
    
  vnl_vector< double > coeffs = Superclass::ComputeCoeffsFromSignal(log_signal, residual);
  
  //int b[5]={0,1,6,15,28};
  //  b[0]=0; b[1]=1; b[2] = 6; b[3]=15; b[4]=28;  l + m + b[ int(double(l)/double(2)) ]
  
  //convert signal RSH coeffs to solid angle odf coeffs
 for( unsigned int i=0; i<OutputPixelType::Dimension; ++i)
  {
    int l = (OutputPixelType::GetLM(i+1))[0];
    coeffs[i] = static_cast<double>( -1 * l * (l+1) * LegendreP( l , 0, 0 ) * coeffs[i]  / ( 8 * vnl_math::pi ) );
  }
    coeffs[0] = 1 / (vcl_sqrt(4 * vnl_math::pi));

/*std::cout<< signal.size() << "   " <<"\n";
for(int i = 0; i < 64; ++i)
{
   std::cout<< "value at : "<< i <<" is:  "<< log_signal.get(i) <<"\n";
}*/
 
return coeffs;
}

template< 
          class TGradientImagePixelType,
          class TOutputPixelType,
          unsigned int TImageDimension,
          class TResidualPercisionType
        >
double SolidAngleOdfReconImageFilter< TGradientImagePixelType, TOutputPixelType, TImageDimension, TResidualPercisionType >
::SmoothingFunction( double x )
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

}//end namespace itk

#endif
