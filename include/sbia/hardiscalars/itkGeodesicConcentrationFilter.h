#ifndef __itkGeodesicConcentrationFilter_h
#define __itkGeodesicConcentrationFilter_h

#include "itkBinaryFunctorImageFilter.h"
#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "itkRSHPixelReorientationOperator.h"
#include "itkMatrix.h"
#include "itkVector.h"
#include "vector.h"

namespace itk
{
  
/** Functor GeodesicConcentrationFilter
 * 
 *  Implements the method described in "Peak Geodesic Concentration: A measure of WM complexity in HARDI", 
 *  which is to appear in the proceedings of MMBIA 2012.
 *
 * 
 * \warning No numeric overflow checking is performed in this filter.
 * 
 **/

namespace Functor{  
  
template< class TInput1, class TInput2, class TOutput, unsigned int ImageDim = 3>
class GeodesicConcentration
{
public:
  
  typedef TInput1                                                  RshPixelType; 
  typedef TInput2                                                  PeaksPixelType;
  typedef TOutput                                                  OutputPixelType;
  typedef typename RshPixelType::ComponentType                     ComponentType;
  
   /** Dimension of the RSH pixel and Peak pixel */
  itkStaticConstMacro(RshDimension, unsigned int, RshPixelType::Dimension );
  itkStaticConstMacro(VecLength, unsigned int, PeaksPixelType::Dimension );
    
  typedef vnl_vector_fixed< ComponentType, ImageDim>               PeakDirectionType;  
  typedef itk::PixelReorientationOperator< RshPixelType >          RotationOperatorType;
  typedef typename   RotationOperatorType::Pointer	           RotationOperatorPointer;

  typedef itk::Matrix< ComponentType, ImageDim, ImageDim>          MatrixType; 
  typedef vnl_vector<ComponentType>                                SamplesVecType;
  
  int                                                              count, nmPeaks;
  ComponentType                                                    epsilon;

GeodesicConcentration() 
{
    count = 0;
    epsilon = .0001;
    m_Rotator = RotationOperatorType::New();
    outVec = 0;
    std::cerr<< " Finished With Constructor" << "\n\n";
}

virtual ~GeodesicConcentration() 
{
 m_Norm = false;
}

inline TOutput operator() ( const RshPixelType & A, const PeaksPixelType & B ) 
{
  
  vnl_vector_fixed< ComponentType , RshDimension>      input;
  RshPixelType C; 
  outVec = 0;
  
 #ifdef VERBOSE 
   std::cout << "Input to GC Operators" << std::endl;
   std::cout << A << std::endl << std::endl;
 #endif
    
   // 'cast' A to a vnl_vector_fixed
   // ...also normalize RSH vector 
  if(m_Norm) 
   {
     C = A;
     C.Normalize();
     for(int i=0; i < RshDimension; ++i) input.put( i , C[i] );
   }
  else
     for(int i=0; i < RshDimension; ++i) input.put( i , A[i] );
     
  if( input.one_norm() >= epsilon )
   {    
       peak.put(0, B[0]);
       peak.put(1, B[1]);
       peak.put(2, B[2]);
       peak = peak.normalize();
       
       if( vcl_fabs( peak.get(2) - 1 ) < epsilon )
         C = A;
       else 
         C = OdfRotation(A, peak);
       
       outVec = vcl_sqrt(vnl_math::pi) * (.2 *C[0] +.255550626*C[3] + .07619047619*C[10] );
   }
   
  if ( outVec > OutputPixelType(0) )
      return static_cast<TOutput>( outVec );
  else 
      return static_cast<TOutput>( 0 );
}

bool operator== (const GeodesicConcentration& other) const
  {
    return !( *this != other );
  }
bool operator!= (const GeodesicConcentration ) const
  {
     return false;
  }

void NormalizationOn(){m_Norm = true;}
void NormalizationOff(){m_Norm = false;}
void GetNormalizationUse(){return m_Norm;}

//Creates a rotation matrix that takes p to q.
inline RshPixelType OdfRotation(RshPixelType odf, PeakDirectionType peak_temp)
{

//Compute rotation matrix
MatrixType A, RotMat;

if( peak_temp.get(2) == -1 )
 {
    RotMat(0,0) = 1; 
    RotMat(0,1) = 0;
    RotMat(0,2) = 0; 
    RotMat(1,0) = 0;
    RotMat(1,1) = 1; 
    RotMat(1,2) = 0; 
    RotMat(2,0) = 0; 
    RotMat(2,1) = 0; 
    RotMat(2,2) = -1;   
    
    return m_Rotator->Reorient(odf, RotMat); 
 }       
    
//Compute Euler Axis;
ComponentType alpha = vcl_acos(peak_temp.get(2));
ComponentType norm_axis = vcl_sin(alpha);
PeakDirectionType euler_axis(-1*peak_temp.get(1)/norm_axis, peak_temp.get(0)/norm_axis, 0);



A(0,0) = euler_axis.get(0); 
A(0,1) = euler_axis.get(1);
A(0,2) = 0; 
A(1,0) = euler_axis.get(1);
A(1,1) = -1*euler_axis.get(0);
A(1,2) = 0; 
A(2,0) = 0; 
A(2,1) = 0; 
A(2,2) = 1;

RotMat(0,0) = 1; 
RotMat(0,1) = 0;
RotMat(0,2) = 0; 
RotMat(1,0) = 0;
RotMat(1,1) = vcl_cos(alpha); 
RotMat(1,2) = -1*norm_axis; 
RotMat(2,0) = 0; 
RotMat(2,1) = norm_axis; 
RotMat(2,2) = vcl_cos(alpha); 

return m_Rotator->Reorient(odf, A*RotMat*A); 
}//end OdfRotation

private:
                                     
PeakDirectionType                               peak;
RotationOperatorPointer                         m_Rotator;
OutputPixelType                                 outVec;
bool m_Norm; //should the RSH functions be normalized. 


}; //end Class

}//end Functor

template <class TInputImage1, class TInputImage2, class TOutputImage> 
class ITK_EXPORT GeodesicConcentrationFilter:
    public
    BinaryFunctorImageFilter< TInputImage1, TInputImage2, TOutputImage,Functor::GeodesicConcentration< typename TInputImage1::PixelType, typename TInputImage2::PixelType, typename TOutputImage::PixelType, TInputImage1::ImageDimension > >
{
public:

  /** Standard class typedefs. */
  typedef GeodesicConcentrationFilter                                    Self;
  typedef BinaryFunctorImageFilter< TInputImage1, TInputImage2, TOutputImage, 
                                    Functor::GeodesicConcentration< typename TInputImage1::PixelType, typename TInputImage2::PixelType, 
                                                         typename TOutputImage::PixelType, TInputImage1::ImageDimension >   >         
                                                                    Superclass;
                           
              
  typedef SmartPointer<Self>                                                  Pointer;
  typedef SmartPointer<const Self>                                            ConstPointer;

  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(GeodesicConcentrationFilter, 
               BinaryFunctorImageFilter);
               
               
void NormalizationOn(){this->GetFunctor().NormalizationOn();}
void NormalizationOff(){this->GetFunctor().NormalizationOff();}
void GetNormalizationUse(){this->GetFunctor().GetNormalizationUse();} 

protected:
  GeodesicConcentrationFilter() {}
  virtual ~GeodesicConcentrationFilter() {}
  
private:

GeodesicConcentrationFilter(const Self&); //purposely not implemented
void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#endif
