#ifndef __itkRshInterpolateImageFunction_txx
#define __itkRshInterpolateImageFunction_txx

#include "itkRshInterpolateImageFunction.h"

/*!
 * \file  itkRSHPixelReorientationOperator.txx
 * \brief -
 *
 * Copyright (c) 2012 University of Pennsylvania. All rights reserved.
 * See COPYING file or https://www.rad.upenn.edu/sbia/software/license.html.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 * 
 *  Implementation of a bspline of order k.
 *  Uses the DeBoor-Cox Algorithm.
 *
 *
 *
 */
namespace itk
{

template< class ImageType, class TComponentType >
TComponentType RshInterpolateFunction::BSplineBasisFunction(const TComponentType x, const int k, const TComponentType begin, const TComponentType end)
    {
       // 0-order BSpline function
        if( k == 0 && begin <= x && x <= end )
            return 1;
        else 
            return 0;
 
        TComponentType alpha = (x - begin) / float(k);
        
        return ((x - begin)/TComponentType(k))*BSplineBasisFunction( x , k-1, begin, end) + ((k+1+begin - x)/TComponentType(k))*BSplineBasisFunction(x, k-1, begin+1, end+1);    
    }
    
template< class ImageType, class TComponentType >    
void RshInterpolateFunction::ComputeBSplineMat(MatrixType &mat, const NeighborhoodIterator a, const NeighborhoodIterator b )
{
    IndexType first, second; 
    
    for( int i=0; i < a.Size() ; ++i )
    {
        first = a.GetIndex(i);
        for( int j = 0; j < b.Size(); ++j) 
        {
            second = first - b.GetIndex(i);
            mat.put(i, j, BSplineBasisFunction(second[0]) * BSplineBasisFunction(second[1]) * BSplineBasisFunction(second[2]) );
        }
    }   
}

template< class ImageType, class TComponentType >
TComponentType RshInterpolateFunction::BSplineFunctionInterpolation(PointType &point, OptimizerType::ParametersType coeff, const NeighborhoodIterator a)
{

   TComponentType accum = 0;
   IndexType first; 
   
   for( int j = 0; j < a.Size(); ++j) 
     {
            first = a.GetIndex(i);
            accum += coeff[j] * BSplineBasisFunction(first[0]) * BSplineBasisFunction(first[1]) * BSplineBasisFunction(first[2]);
     }

   return accum;
}

    
} // end namespace itk

#endif
