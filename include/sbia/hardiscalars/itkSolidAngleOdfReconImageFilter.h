#ifndef __itkSolidAngleOdfReconImageFilter_h_
#define __itkSolidAngleOdfReconImageFilter_h_

#include "itkSymRealSphericalHarmonicRep.h"
#include "itkSymRshReconImageFilter.h"

namespace itk
{

template< 
          class TGradientImagePixelType,
          class TOutputPixelType = SymRealSphericalHarmonicRep< float, 4 >,
          unsigned int TImageDimension=3,
          class TResidualPercisionType = float
        >
class ITK_EXPORT SolidAngleOdfReconImageFilter :
  public SymRshReconImageFilter<  TGradientImagePixelType,
                                  TOutputPixelType, TImageDimension, TResidualPercisionType >

{
public:

  typedef SolidAngleOdfReconImageFilter           Self;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;

  typedef SymRshReconImageFilter< TGradientImagePixelType,
                                  TOutputPixelType, TImageDimension, TResidualPercisionType >
                                                  Superclass;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(SolidAngleOdfReconImageFilter, SymRshReconImageFilter);

  /** pass through typdefs from superclass. */
  typedef typename Superclass::GradientPixelType                GradientPixelType;
  typedef typename Superclass::ResidualPixelType                ResidualPixelType;

  typedef typename Superclass::OutputImageType                  OutputImageType;
  typedef typename Superclass::OutputPixelType                  OutputPixelType;
  typedef typename Superclass::OutputImageRegionType            OutputImageRegionType;
  typedef typename Superclass::GradientImagesType               GradientImagesType;
  typedef typename Superclass::GradientDirectionContainerType   GradientDirectionContainerType;
  
    
    /** Set/Get amount of 'smoothing' applied to signal  */
    itkSetMacro(Delta, double );
    itkGetMacro(Delta, double );
    
protected:
    SolidAngleOdfReconImageFilter();
  ~SolidAngleOdfReconImageFilter() {};
  

  /** just overload the computeCoeffsFromSignal method. */
  virtual vnl_vector< double >  ComputeCoeffsFromSignal( vnl_vector< double >,ResidualPixelType& );
  
private:

//The parameter that determines how much "smoothing" applied to the signal
double                                                          m_Delta;

/** function that smooths DWI signal           **/
double SmoothingFunction( double );

};

} //end itkNamespace

#include "itkSolidAngleOdfReconImageFilter.txx"

#endif
