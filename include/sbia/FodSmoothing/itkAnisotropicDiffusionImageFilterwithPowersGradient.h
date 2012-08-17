/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAnisotropicDiffusionImageFilterwithPowersGradient.h,v $
  Language:  C++
  Date:      $Date: 2012-06-25 14:31:35 -0400 (Mon, 25 Jun 2012) $
  Version:   $Revision: 38 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAnisotropicDiffusionImageFilterwithPowersGradient_h
#define __itkAnisotropicDiffusionImageFilterwithPowersGradient_h

#include "itkExceptionObject.h"
#include "itkAnisotropicDiffusionImageFilter.h"
#include "itkAnisotropicDiffusionImageFunctionwithPowersGradient.h"

namespace itk {

/** \class itkAnisotropicDiffusionImageFilterwithPowersGradient
 *
 * This filter performs anisotropic diffusion on a vector itk::Image using the
 * anisotropic diffusion function implemented implemented in
 * itkVectorGradientNDAnisotropicDiffusionFunction.  For detailed information on
 * anisotropic diffusion see itkAnisotropicDiffusionFunction,
 * itkVectorGradientNDAnisotropicDiffusionFunction, and
 * itkGradientAnisotropicDiffusionFunction.
 * 
 * \par Inputs and Outputs
 * The input to this filter must be an itk::Image with pixel
 * type which is either an itk::Vector, or a subclass of an itk::Vector.
 * Additionally, the component type of the vector should be a numerical type
 * (float or double, or a user defined type which correctly defines
 * arithmetic operations with floating point accuracy).  The output image type
 * also has these requirements.
 *
 * \par Parameters
 * Please read all the documentation found in
 * AnisotropicDiffusionImageFilter and AnisotropicDiffusionFunction.  Also see
 * VectorGradientNDAnisotropicDiffusionFunction.
 *
 * The maximum allowable time step for this filter is 1/2^N, where N is the
 * dimensionality of the image.  For 2D images any value below 0.250 is stable,
 * and for 3D images, any value below 0.125 is stable.
 *
 * \ingroup ImageEnhancement
 */
template <class TInputImage, class TOutputImage>
class ITK_EXPORT AnisotropicDiffusionImageFilterwithPowersGradient
  : public AnisotropicDiffusionImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef AnisotropicDiffusionImageFilterwithPowersGradient Self;
  typedef AnisotropicDiffusionImageFilter<TInputImage, TOutputImage>
                                                        Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Instantiation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information. */
  itkTypeMacro(AnisotropicDiffusionImageFilterwithPowersGradient,
               AnisotropicDiffusionImageFilter);
  
  /** Extract information from the superclass. */
  typedef typename Superclass::UpdateBufferType UpdateBufferType;
  
   /** Determine the image dimension from the  superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      Superclass::ImageDimension );
		      
  typedef typename TInputImage::InternalPixelType   PixelType;
  
  /** Type of a value in a vector (double, float, etc.) */
  typedef typename PixelType::ValueType  ScalarValueType;
  
/* const unsigned int drivingComponents	;      */
  virtual void SetDrivingComponents(const unsigned int num)
    {
    typename AnisotropicDiffusionImageFunctionwithPowersGradient<UpdateBufferType>::Pointer p
      = AnisotropicDiffusionImageFunctionwithPowersGradient<UpdateBufferType>::New();
    this->SetDifferenceFunction(p);
    p->SetDrivingComponents(num);
    std::cout<< "driving element is" << num << std::endl;
    }
    
  
#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<typename TInputImage::PixelType::ValueType>));
  itkConceptMacro(OutputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<typename TOutputImage::PixelType::ValueType>));
  /** End concept checking */
#endif

protected:
  AnisotropicDiffusionImageFilterwithPowersGradient()
  {
  }
  
  ~AnisotropicDiffusionImageFilterwithPowersGradient() {}

private:
  AnisotropicDiffusionImageFilterwithPowersGradient(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namspace itk

#endif
