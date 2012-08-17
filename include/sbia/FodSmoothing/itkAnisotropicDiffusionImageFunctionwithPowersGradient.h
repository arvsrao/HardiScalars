/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAnisotropicDiffusionImageFunctionwithPowersGradient.h,v $
  Language:  C++
  Date:      $Date: 2012-06-25 14:31:35 -0400 (Mon, 25 Jun 2012) $
  Version:   $Revision: 38 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAnisotropicDiffusionImageFunctionwithPowersGradient_h
#define __itkAnisotropicDiffusionImageFunctionwithPowersGradient_h

#include "itkVectorAnisotropicDiffusionFunction.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkVectorNeighborhoodInnerProduct.h"
#include "itkDerivativeOperator.h"
#include "itkGaussianOperator.h"
#include "vnl/vnl_math.h"
#include <math.h>


namespace itk {

/** \class AnisotropicDiffusionImageFunctionwithPowersGradient
 *
 * This class is a simple extension of the
 * GradientNDAnisotropicDiffusionFunction to pixel types of multiple
 * components.  Vector components are diffused separately, but diffusion of
 * each component is limited by a conductance term which depends on all
 * components.
 *
 * For more information, please see GradientNDAnisotropicDiffusionFunction.
 *
 * \sa GradientNDAnisotropicDiffusionFunction
 * \sa VectorCurvatureNDAnisotropicDiffusionFunction
 * \sa AnisotropicDiffusionFunction
 */ 
template <class TImage>
class ITK_EXPORT AnisotropicDiffusionImageFunctionwithPowersGradient :
    public VectorAnisotropicDiffusionFunction<TImage>
{
public:
  /** Standard class typedefs. */
  typedef AnisotropicDiffusionImageFunctionwithPowersGradient Self;
  typedef VectorAnisotropicDiffusionFunction<TImage>   Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( AnisotropicDiffusionImageFunctionwithPowersGradient,
                ScalarAnisotropicDiffusionFunction );
  
  /** Inherit some parameters from the superclass type. */
  typedef typename Superclass::ImageType        ImageType;
  typedef typename Superclass::PixelType        PixelType;
  typedef typename Superclass::TimeStepType     TimeStepType;
  typedef typename Superclass::RadiusType       RadiusType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  typedef typename Superclass::FloatOffsetType  FloatOffsetType;

  /** Extract vector and image dimension from superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      Superclass::ImageDimension );
  itkStaticConstMacro(VectorDimension, unsigned int,
                      Superclass::VectorDimension );


  /** Type of a value in a vector (double, float, etc.) */
  typedef typename PixelType::ValueType  ScalarValueType;
  
  
    /** Compute the equation value. */
  virtual PixelType ComputeUpdate(const NeighborhoodType &neighborhood,
                                  void * globalData,
                                  const FloatOffsetType& offset = FloatOffsetType(0.0)
    );


  /** This method is called prior to each iteration of the solver. */
  virtual void InitializeIteration()
    {
    m_K = this->GetAverageGradientMagnitudeSquared() * this->GetConductanceParameter() *
      this->GetConductanceParameter() * -2.0f;
    }

  unsigned int drivingComponents	;      
  virtual void SetDrivingComponents(const unsigned int num)
    {
     drivingComponents=num;
    }
    
   
  virtual unsigned int GetDrivingComponents()
    {
     return drivingComponents;
    }    


const unsigned int max (PixelType a) {
PixelType b;
b.Fill(0);
b[0]=a[5];
for (unsigned int i=1;i<6;i++)
{
b[1]+=vnl_math_sqr(a[drivingComponents+i]);
}
b[1]=sqrt(b[1]);
for (unsigned int i=6;i<15;i++)
{
b[2]+=vnl_math_sqr(a[drivingComponents+i]);
}
b[2]=sqrt(b[2]);
for (unsigned int i=15;i<28;i++)
{
b[3]+=vnl_math_sqr(a[drivingComponents+i]);
}
b[3]=sqrt(b[3]);
for (unsigned int i=28;i<45;i++)
{
b[4]+=vnl_math_sqr(a[drivingComponents+i]);
}
b[4]=sqrt(b[4]);

double temp = b.GetElement(0);
unsigned int index=0;
for(unsigned int i=1;i<this->GetDrivingComponents();i++)
{
if(temp<b.GetElement(i))
{
  temp = b.GetElement(i); 
  index=  i; 
  }
  }
  
  return  index;
} 

   
      
protected:
  AnisotropicDiffusionImageFunctionwithPowersGradient();
  ~AnisotropicDiffusionImageFunctionwithPowersGradient() {}

private:
  AnisotropicDiffusionImageFunctionwithPowersGradient(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Inner product function. */
  VectorNeighborhoodInnerProduct<ImageType> m_InnerProduct;

  /** Slices for the ND neighborhood. */
  std::slice  x_slice[ImageDimension];
  std::slice xa_slice[ImageDimension][ImageDimension];
  std::slice xd_slice[ImageDimension][ImageDimension];

  /** Derivative operators. */
  DerivativeOperator<ScalarValueType,
                     itkGetStaticConstMacro(ImageDimension)> dx_op;
		     
 GaussianOperator<ScalarValueType,
                     itkGetStaticConstMacro(ImageDimension)> g_op;
		     

  /** Modified global average gradient magnitude term. */
  ScalarValueType m_K;
  
  static double m_MIN_NORM;
  
  unsigned long int m_Stride[ImageDimension];
  unsigned long int m_Center;
  
};


  
}// end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAnisotropicDiffusionImageFunctionwithPowersGradient.txx"
#endif

#endif
