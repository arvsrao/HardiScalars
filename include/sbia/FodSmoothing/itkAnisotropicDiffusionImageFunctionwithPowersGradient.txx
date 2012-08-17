/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAnisotropicDiffusionImageFunctionwithPowersGradient.txx,v $
  Language:  C++
  Date:      $Date: 2012-06-25 14:31:35 -0400 (Mon, 25 Jun 2012) $
  Version:   $Revision: 38 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkAnisotropicDiffusionImageFunctionwithPowersGradient_txx
#define __itkAnisotropicDiffusionImageFunctionwithPowersGradient_txx


namespace itk {

template<class TImage>
double AnisotropicDiffusionImageFunctionwithPowersGradient<TImage>
::m_MIN_NORM = 1.0e-10;
  
template<class TImage>
AnisotropicDiffusionImageFunctionwithPowersGradient<TImage>
::AnisotropicDiffusionImageFunctionwithPowersGradient()
{
  unsigned int i, j;
  RadiusType r;

  for (i = 0; i < ImageDimension; ++i)
    {
    r[i] = 1;
    }
  this->SetRadius(r);

  // Dummy neighborhood used to set up the slices.
  Neighborhood<PixelType, ImageDimension> it;
  it.SetRadius(r);
  
  // Slice the neighborhood
  m_Center =  it.Size() / 2;

  for (i = 0; i< ImageDimension; ++i)
    { m_Stride[i] = it.GetStride(i); }

  for (i = 0; i< ImageDimension; ++i)
    { x_slice[i]  = std::slice( m_Center - m_Stride[i], 3, m_Stride[i]); }
  
  for (i = 0; i< ImageDimension; ++i)
    {
    for (j = 0; j < ImageDimension; ++j)
      {
      // For taking derivatives in the i direction that are offset one
      // pixel in the j direction.
      xa_slice[i][j]
        = std::slice((m_Center + m_Stride[j])-m_Stride[i], 3, m_Stride[i]); 
      xd_slice[i][j]
        = std::slice((m_Center - m_Stride[j])-m_Stride[i], 3, m_Stride[i]);
      }
    }
  
  // Allocate the derivative operator.
  dx_op.SetDirection(0); // Not relelevant, we'll apply in a slice-based
                         // fashion 
  dx_op.SetOrder(1);
  dx_op.CreateDirectional();
  g_op.SetVariance(1);
  g_op.SetDirection(0);
  g_op.CreateDirectional();
  g_op.SetRadius(it.GetRadius());
  
  
  
}

template<class TImage>
typename AnisotropicDiffusionImageFunctionwithPowersGradient<TImage>::PixelType
AnisotropicDiffusionImageFunctionwithPowersGradient<TImage>
::ComputeUpdate(const NeighborhoodType &it, void *,
                const FloatOffsetType&)
{

  unsigned int i, j;

  double accum;
  double accum_d;
  double Cx;
  double Cxd;
  
  // ValueType of each of the components of PixelType - scalar 
  PixelType delta;
  PixelType dx_forward;
  PixelType dx_backward;
  PixelType dx[ImageDimension];
  ScalarValueType dx_aug;
  ScalarValueType dx_dim;
  
  delta = NumericTraits<ScalarValueType>::Zero;
     
  unsigned int chosen_index=max(it.GetPixel(m_Center));
  
  
  
   // Calculate the centralized derivatives for each dimension.
  for (i = 0; i < ImageDimension; i++)
    {
    dx[i]  =  (it.GetPixel(m_Center + m_Stride[i])-it.GetPixel(m_Center - m_Stride[i]))/2.0f;
    dx[i] *= this->m_ScaleCoefficients[i];
    }


  for (i = 0; i < ImageDimension; i++)
    {
    // ``Half'' directional derivatives
    dx_forward = it.GetPixel(m_Center + m_Stride[i])- it.GetPixel(m_Center);
    dx_forward *= this->m_ScaleCoefficients[i];
    dx_backward =  it.GetPixel(m_Center)  - it.GetPixel(m_Center - m_Stride[i]);
    dx_backward *= this->m_ScaleCoefficients[i];

    // Calculate the conductance terms.  Conductance varies with each
    // dimension because the gradient magnitude approximation is different
    // along each  dimension.
    
    
    
    accum   = 0.0;
    accum_d = 0.0;
    for (j = 0; j < ImageDimension; j++)
      {
      if (j != i)
        {
	
	
        dx_aug = (it.GetPixel(m_Center + m_Stride[i] + m_Stride[j]).GetElement(chosen_index) -
                  it.GetPixel(m_Center + m_Stride[i] - m_Stride[j]).GetElement(chosen_index) ) / 2.0f;
        dx_aug *= this->m_ScaleCoefficients[j];
        dx_dim = (it.GetPixel(m_Center - m_Stride[i] + m_Stride[j]).GetElement(chosen_index) -
                  it.GetPixel(m_Center - m_Stride[i] - m_Stride[j]).GetElement(chosen_index) ) /2.0f;
        dx_dim *= this->m_ScaleCoefficients[j];
        accum += 0.25f * vnl_math_sqr( dx[j].GetElement(chosen_index) + dx_aug );
        accum_d += 0.25f * vnl_math_sqr( dx[j].GetElement(chosen_index) + dx_dim );
        }
      }
      
    if (m_K == 0.0)
      {
      Cx = 0.0;
      Cxd = 0.0;
      }
    else
      {
      Cx = vcl_exp(( vnl_math_sqr( dx_forward.GetElement(chosen_index) ) + accum)  / m_K );
      Cxd= vcl_exp(( vnl_math_sqr( dx_backward.GetElement(chosen_index)) + accum_d)/ m_K );
      }

    // Conductance modified first order derivatives.
    dx_forward  = dx_forward * Cx;
    dx_backward = dx_backward * Cxd;

    // Conductance modified second order derivative.
    delta += dx_forward - dx_backward;
    }
    
    
   if(it.GetIndex().GetElement(0)==50 && it.GetIndex().GetElement(1)==50 && it.GetIndex().GetElement(2)==50)
  {
  std::cout<< "driving element is" << chosen_index << std::endl;//"," << it.GetPixel(50) << "," << it.GetPixel(50).GetElement(this->GetDrivingComponent()) << std::endl;
  std::cout<< "delta update value is" << delta  << std::endl;
  } 
  
  
  return delta ;
  
}

} // end namespace itk

#endif
