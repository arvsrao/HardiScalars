/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScalarGradientNDAnisotropicVectorDiffusionFunction..txx,v $
  Language:  C++
  Date:      $Date: 2012-06-25 14:31:35 -0400 (Mon, 25 Jun 2012) $
  Version:   $Revision: 38 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkScalarGradientNDAnisotropicVectorDiffusionFunction_txx
#define __itkScalarGradientNDAnisotropicVectorDiffusionFunction_txx

namespace itk {

template<class TImage>
double ScalarGradientNDAnisotropicVectorDiffusionFunction<TImage>
::m_MIN_NORM = 1.0e-10;
  
template<class TImage>
ScalarGradientNDAnisotropicVectorDiffusionFunction<TImage>
::ScalarGradientNDAnisotropicVectorDiffusionFunction()
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
}

template<class TImage>
typename ScalarGradientNDAnisotropicVectorDiffusionFunction<TImage>::PixelType
ScalarGradientNDAnisotropicVectorDiffusionFunction<TImage>
::ComputeUpdate(const NeighborhoodType &it, void *,
                const FloatOffsetType&)
{

  unsigned int i, j;

  double accum;
  double accum_d;
  double Cx;
  double Cxd;
  
  // ValueType of each of the components of PixelType - scalar 
  ScalarValueType delta;
  ScalarValueType dx_forward;
  ScalarValueType dx_backward;
  ScalarValueType dx[ImageDimension];
  ScalarValueType dx_aug;
  ScalarValueType dx_dim;
  PixelType delta_return;

  delta = NumericTraits<ScalarValueType>::Zero;
  
  // Calculate the centralized derivatives for each dimension.
  for (i = 0; i < ImageDimension; i++)
    {
    dx[i]  =  (it.GetPixel(m_Center + m_Stride[i]).GetElement(this->GetDrivingComponent())-it.GetPixel(m_Center - m_Stride[i]).GetElement(this->GetDrivingComponent()))/2.0f;
    dx[i] *= this->m_ScaleCoefficients[i];
    }

  for (i = 0; i < ImageDimension; i++)
    {
    // ``Half'' directional derivatives
    dx_forward = it.GetPixel(m_Center + m_Stride[i]).GetElement(this->GetDrivingComponent())
      - it.GetPixel(m_Center).GetElement(this->GetDrivingComponent());
    dx_forward *= this->m_ScaleCoefficients[i];
    dx_backward =  it.GetPixel(m_Center).GetElement(this->GetDrivingComponent())
      - it.GetPixel(m_Center - m_Stride[i]).GetElement(this->GetDrivingComponent());
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
        dx_aug = (it.GetPixel(m_Center + m_Stride[i] + m_Stride[j]).GetElement(this->GetDrivingComponent()) -
                  it.GetPixel(m_Center + m_Stride[i] - m_Stride[j]).GetElement(this->GetDrivingComponent()) ) / 2.0f;
        dx_aug *= this->m_ScaleCoefficients[j];
        dx_dim = (it.GetPixel(m_Center - m_Stride[i] + m_Stride[j]).GetElement(this->GetDrivingComponent()) -
                  it.GetPixel(m_Center - m_Stride[i] - m_Stride[j]).GetElement(this->GetDrivingComponent()) ) /2.0f;
        dx_dim *= this->m_ScaleCoefficients[j];
        accum += 0.25f * vnl_math_sqr( dx[j] + dx_aug );
        accum_d += 0.25f * vnl_math_sqr( dx[j] + dx_dim );
        }
      }
      
    if (m_K == 0.0)
      {
      Cx = 0.0;
      Cxd = 0.0;
      }
    else
      {
      Cx = vcl_exp(( vnl_math_sqr( dx_forward ) + accum)  / m_K );
      Cxd= vcl_exp(( vnl_math_sqr( dx_backward) + accum_d)/ m_K );
      }

    // Conductance modified first order derivatives.
    dx_forward  = dx_forward * Cx;
    dx_backward = dx_backward * Cxd;

    // Conductance modified second order derivative.
    delta += dx_forward - dx_backward;
    }
  for(unsigned int k =0; k < VectorDimension; k++)
  {
	  delta_return[k]=delta;
  }
/*   if(it.GetIndex().GetElement(0)==50 && it.GetIndex().GetElement(1)==50 && it.GetIndex().GetElement(2)==50)
  {
  std::cout<< "driving element is" << drivingComponent << std::endl;//"," << it.GetPixel(50) << "," << it.GetPixel(50).GetElement(this->GetDrivingComponent()) << std::endl;
  std::cout<< "delta update value is" << delta_return << "," << delta << std::endl;
  } */
  return delta_return ;
  
  ///////////////////////////////////
  /* unsigned int i, j, k;
  PixelType delta;

  double GradMag;
  double GradMag_d;
  double Cx[ImageDimension];
  double Cxd[ImageDimension];

  // Remember: PixelType is a Vector of length VectorDimension.
  PixelType dx_forward[ImageDimension];
  PixelType dx_backward[ImageDimension];
  PixelType dx[ImageDimension];
  PixelType dx_aug;
  PixelType dx_dim;

  // Calculate the directional and centralized derivatives.
  for (i = 0; i < ImageDimension; i++)
    {
    // ``Half'' derivatives
    dx_forward[i] = it.GetPixel(m_Center + m_Stride[i])
      - it.GetPixel(m_Center);
    dx_forward[i] = dx_forward[i]  * this->m_ScaleCoefficients[i];
    dx_backward[i]=  it.GetPixel(m_Center)
      - it.GetPixel(m_Center - m_Stride[i]);
    dx_backward[i] = dx_backward[i] * this->m_ScaleCoefficients[i];

    // Centralized differences
    dx[i]      = m_InnerProduct(x_slice[i], it, dx_op);
    dx[i] = dx[i] * this->m_ScaleCoefficients[i];
    }
    
    
  // Calculate the conductance term for each dimension.
  for (i = 0; i < ImageDimension; i++)
    {
    // Calculate gradient magnitude approximation in this
    // dimension linked (summed) across the vector components.
    GradMag   = 0.0;
    GradMag_d = 0.0;
    
    for (k =0; k < VectorDimension; k++) 
      {
      GradMag += vnl_math_sqr( dx_forward[i][k] );
      GradMag_d += vnl_math_sqr( dx_backward[i][k] );

      for (j = 0; j < ImageDimension; j++)
        {
        if ( j != i)
          {
          dx_aug  = m_InnerProduct(xa_slice[j][i], it, dx_op);
          dx_aug = dx_aug * this->m_ScaleCoefficients[j];
          dx_dim  = m_InnerProduct(xd_slice[j][i], it, dx_op);
          dx_dim = dx_dim * this->m_ScaleCoefficients[j];
          GradMag += 0.25f * vnl_math_sqr( dx[j][k]+dx_aug[k] );
          GradMag_d += 0.25f * vnl_math_sqr( dx[j][k]+dx_dim[k] );
          }
        }
      }
      
    if (m_K == 0.0)
      {
      Cx[i] = 0.0;
      Cxd[i] = 0.0;
      }
    else
      {
      Cx[i]  = vcl_exp( GradMag   / m_K );
      Cxd[i] = vcl_exp( GradMag_d / m_K ); 
      }
    }

  // Compute update value  
  for (k = 0; k < VectorDimension; k++)
    {
    delta[k] = NumericTraits<ScalarValueType>::Zero;
      
    for (i = 0; i < ImageDimension; ++i)
      {
      dx_forward[i][k] *= Cx[i];
      dx_backward[i][k] *= Cxd[i];
      delta[k] += dx_forward[i][k] - dx_backward[i][k];
      }
    }
      
  return delta; */
}

} // end namespace itk

#endif
