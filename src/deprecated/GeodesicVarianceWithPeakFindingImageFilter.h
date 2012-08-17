#ifndef __itkGeodesicVarianceWithPeakFindingImageFilter_h
#define __itkGeodesicVarianceWithPeakFindingImageFilter_h

#include "itkBinaryFunctorImageFilter.h"
#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "itkRSHPixelReorientationOperator.h"
#include "itkPeakFindingCalculator.h"
#include "itkMatrix.h"
#include "itkVector.h"
#include "vector.h"
#include "sphere1000points.h"


namespace itk
{
  
/** \class GeodesicVarianceWithPeakFindingImageFilter
 * #include <itkMatrix.h> * Implements an operator for calculating GeodesicVarianceWithPeakFinding scalar measures on
 * ODF images. 
 *
 * 
 * \warning No numeric overflow checking is performed in this filter.
 * 
 **/

namespace Functor{  
  
template< unsigned int ImageDim = 3, class TInput1, class TInput2, class TOutput>
class GeodesicVarianceWithPeakFinding
{
public:
  
  typedef TInput1                                                     MaskPixelType;
  typedef TInput2                                                     OdfPixelType; 
  typedef TOutput                                                     OutputPixelType;
  
  typedef typename OdfPixelType::ComponentType                        ComponentType;
  
  typedef vnl_vector_fixed<ComponentType, ImageDim>                   PeakDirectionType;  
  typedef itk::PixelReorientationOperator< OdfPixelType >             RotationOperatorType;
  
  typedef itk::PeakFindingCalculator< OdfPixelType, ComponentType , ComponentType >  
                                                                      PeakFindingCalculatorType;
  
  typedef typename   RotationOperatorType::Pointer	              RotationOperatorPointer;
  typedef typename   PeakFindingCalculatorType::Pointer               PeakFindingCalculatorPointer;
  
  typedef itk::Matrix< ComponentType, ImageDim, ImageDim>             MatrixType; 
  
  itkStaticConstMacro( MaxOrder,  unsigned int, OdfPixelType::MaxOrder  );
  itkStaticConstMacro( Dimension, unsigned int, OdfPixelType::Dimension );
 
  typedef typename PeakFindingCalculatorType::InputType               PeakFindingVectorType;
   
  typedef vnl_matrix< ComponentType >                                          VnlMatrixType;
  typedef vnl_vector_fixed< ComponentType, 10000>                              VectorType;
  

GeodesicVarianceWithPeakFinding() 
{

    m_Rotator = RotationOperatorType::New();
    m_PeakFinder = PeakFindingCalculatorType::New();
    mat.set_size(10000, Dimension);
    limit[0]=5; limit[1]=5; limit[2]=5; limit[3]=5;
    geoVar.SetSize( 10 ); isotropic.SetSize(10);
    for(unsigned int k = 0; k < 10; ++k){ isotropic[k]= 0; geoVar[k]=0; }
 
    //Fill theta and phi vectors
    for(int i=0; i < 10000; ++i)
    {
         theta.put( i, vcl_acos( points10000Z[i] ) );

         if( (points10000Z[i] != 1) & (points10000Z[i] != -1) )
         {
            if( (points10000X[i] == 0) & (points10000Y[i] < 0) )
  	        phi.put(i, -1 * vnl_math::pi);
            else if( (points10000X[i]==0) & (points10000Y[i] > 0) )
                phi.put(i, vnl_math::pi);
            else if ( (points10000X[i] < 0) & (points10000Y[i] > 0) )
                phi.put(i, vnl_math::pi/2 + vcl_atan( points10000Y[i] / points10000X[i]) ); 
	    else if ( (points10000X[i] < 0) & (points10000Y[i] < 0) )
                phi.put(i, vnl_math::pi + vcl_atan( points10000Y[i] / points10000X[i] ) );
            else
                phi.put(i, vcl_atan( points10000Y[i] / points10000X[i]) );
         }
         else 
             phi.put(i, 0);
    }

   //create mat matrix of sampled SH basis functions
     for(int i=0; i < 10000; ++i)
     {
        for(int j = 1; j <= Dimension ; ++j)
        { 
            mat.put(i, j-1, TInput2::Y(j, theta.get(i), phi.get(i) ) );
        }       
     }
}

virtual ~GeodesicVarianceWithPeakFinding() {}


inline TOutput operator() ( const TInput1 & A, const TInput2 & B ) 
{

//double maxMin[2]={1, 0.1};

//fill shCoefs and MinShCoefs vectors
for(int i = 0; i < Dimension; ++i)
 {
   shCoefs[i] = B[i];
   if(i==0)
      MinShCoefs[0] = .2 - B[0];
   else
      MinShCoefs[i] = -1*B[i];
 }
   
if( A!=1 || shCoefs.GetNorm()==0)
 {
  
   isotropic[0] = GeneralFractionalAnisotropyFunction(B);
   isotropic[1] = 1 - GeoVarianceFunction( B , vnl_math::pi / (double) 2 );
   
   for(unsigned int k = 2; k <= 4; ++k){ isotropic[k]=isotropic[1]; }
   
   isotropic[5] = isotropic[1] / isotropic[0];         //  GV_restricted / GFA 
   isotropic[6] = 1 / isotropic[5];                    //  recipricol
   isotropic[7] = ( isotropic[0] + isotropic[1] ) / double(2);         //  average of GFA and unrestricted GV for lobe #1
   isotropic[8] = isotropic[7];
   isotropic[9] = isotropic[1];                        //  average of restricted GV
   
   
   return static_cast<TOutput>( isotropic );
 }
else
   { 
      odPeaks.clear(); odMins.clear();
      nmPeaks = 0; nmMins = 0;
      limit[0]=5; limit[1]=5;
      
      m_PeakFinder->SetAcceleration( 0 );
      m_PeakFinder->SetPeakFiltering( 1 );
      //m_PeakFinder->SetRangeOfInterest( maxMin );
      
      //Get Peaks
      m_PeakFinder->SetCoefficients( shCoefs );
      nmPeaks = m_PeakFinder->GetRawPeaks( odPeaks );
      
      //Get Mins
      m_PeakFinder->SetCoefficients( MinShCoefs );
      nmMins = m_PeakFinder->GetRawPeaks( odMins ); 
        
     //Threshold to determine the number of 'genuine' peaks and calculate GV   
     for(int i = 0; i < int(odPeaks.size()); ++i)
      { 
          if ( i % 4 == 3 & odPeaks[i] >= .085 ) 
           {   
              ++nmPeaks;
           }
      }
       
      //Threshold to determine the number of 'genuine' Mins and calculate GV   
     for(int i = 0; i < int(odMins.size()); ++i)
      { 
          if ( i % 4 == 3 & odMins[i] >= .92 ) 
           {   
              ++nmMins;
           }
      }
        
      //Calculate Discrete GFA  
      geoVar[0] = GeneralFractionalAnisotropyFunction(B); //geoVar[0] = nmPeaks; 
        
     //determine limits of integration (i.e. boundaries of lobes)
     for (int i = 1; i <= vnl_math_min(int(2), nmPeaks); ++i )
      {
          peak.put( 0, double(odPeaks[4*(i-1)]) );
          peak.put( 1, double(odPeaks[1+4*(i-1)]) ); 
          peak.put( 2, double(odPeaks[2+4*(i-1)]) );
                 
          for(int j = 1; j <= nmMins; ++j)
           {
             min.put( 0, double(odMins[4*(j-1)]) );
             min.put( 1, double(odMins[1 + 4*(j-1)]) );
             min.put( 2, double(odMins[2 + 4*(j-1)]) );
             
             if( limit[i-1] > vcl_acos(dot_product(peak, min)) )      
                 limit[i-1] = vcl_acos(dot_product(peak, min));     
           }
      }
      
     //std::cout<< limit[ <<"\n";
     //limit the number of discoverable fibers to 2. 
     for( int j = 1; j <=vnl_math_min(int(2), nmPeaks); ++j )
      { 
         
           peak.put( 0, double(odPeaks[4*(j-1)]) );
           peak.put( 1, double(odPeaks[1+4*(j-1)]) ); 
           peak.put( 2, double(odPeaks[2+4*(j-1)]) );
         
         if ( nmPeaks < 2 )
          {
             geoVar[1] = 1 - GeoVarianceFunction( OdfRotation(B, peak), vnl_math::pi / (double) 2 );
             geoVar[3] = geoVar[1];
             geoVar[2] = 0; geoVar[4] = 0; 
             break;
          }    
                        
           geoVar[j] = 1 - GeoVarianceFunction( OdfRotation(B, peak), limit[j-1] );
           geoVar[j+2] = 1 - GeoVarianceFunction( OdfRotation(B, peak), vnl_math::pi / (double) 2 );
      }
    
      /*//Fill in zeros for non-peak indices up to j=4, leaving j=5,6 open for mixtures of computed GV's
      for( int j = nmPeaks + 1; j <=2; ++j )
       {
          geoVar[j] = 0; 
       }
     */
      //Compute mixtures of GeoVariances
          geoVar[5] = geoVar[1] / geoVar[0];                  //  GV_restricted / GFA 
          geoVar[6] = 1 / geoVar[5];                          //  recipricol
          geoVar[7] = (geoVar[0] + geoVar[3])/double(2);      //  average of GFA and unrestricted GV for lobe #1
          geoVar[8] = (geoVar[0] + geoVar[1])/double(2);      //  average of GFA and restricted GV for lobe #1
          geoVar[9] = (geoVar[1] + geoVar[2])/double(2);      //  average of restricted GV
         
      return static_cast<TOutput>( geoVar );
   }
}

bool operator== (const GeodesicVarianceWithPeakFinding& other) const
  {
    return !( *this != other );
  }
bool operator!= (const GeodesicVarianceWithPeakFinding ) const
  {
     return false;
  }

private:
/*Private Member functions*/

int                                                       nmPeaks, nmMins;
std::vector<float>                                        odPeaks, odMins;
double                                                    limit[4];
PeakFindingVectorType                                     shCoefs, MinShCoefs;
PeakDirectionType                                         peak, min;
OutputPixelType                                           geoVar, isotropic;

VnlMatrixType mat; 
VectorType theta; 
VectorType phi;
vnl_matrix_fixed<double, 10000, 3> points;

RotationOperatorPointer                                             m_Rotator;
PeakFindingCalculatorPointer                                        m_PeakFinder;
      
float GeoVarianceFunction(const TInput2 & C, double limit)
  {
     vnl_vector_fixed<double , Dimension>      input;
     double                                    OneNorm=0;
     double                                    sum=0;
     
     //Cast C to a vnl_vector_fixed
     for(int i=0; i < Dimension ; ++i)
      {
        input.put( i , C[i] );
      }
      
     VectorType                                samples = mat * input;
                   
     for(int i=0; i <10000; ++i)
      {
        if( theta.get(i) <= limit )
           OneNorm += samples.get(i) * vcl_sin( theta.get(i) );
      }

     for(int j = 0; j < 10000; ++j)
      {
	 if( theta.get(j) <= limit )
            // sum += samples.get(j) * vcl_pow( vcl_cos( 2 * theta.get(j) /(double)4 ), (double)4 ) * vcl_sin(theta.get(j));
	     sum += samples.get(j) * theta.get(j) *theta.get(j) * vcl_sin(theta.get(j));
      }
      
      return ( sum/(OneNorm * limit * limit ) );
  }

float GeneralFractionalAnisotropyFunction( const TInput2 & A )
{  
 
   vnl_vector_fixed<double , Dimension> input;
   for(int i=0; i < Dimension ; ++i)
    {
     input.put( i , A[i] );
    }
   
   VectorType                          samples = mat * input; 
   double                              sstd, rms, mean;
   
   mean = samples.mean();
   rms  = dot_product(samples, samples);
   sstd = rms - 2*mean*samples.sum() + mean * mean * double(10000); 

   sstd /= double(9999);
   rms  /= double(10000);

   // std::cout<< "GFA value is: " << vcl_sqrt( sstd / rms ) <<"\n";
        
   return vcl_sqrt( sstd / rms );
}
 
//Creates a rotation matrix that takes p to q.
OdfPixelType OdfRotation(OdfPixelType odf, PeakDirectionType p)
{

float alpha = vcl_acos(p.get(2));

//Compute rotation matrix
MatrixType R_x, R_y;
R_x(0,0) = 1; 
R_x(0,1) = 0;
R_x(0,2) = 0; 
R_x(1,0) = 0;
R_x(1,1) =  vcl_cos(alpha); 
R_x(1,2) = -vcl_sin(alpha); 
R_x(2,0) = 0; 
R_x(2,1) = vcl_sin(alpha); 
R_x(2,2) = vcl_cos(alpha); 

R_y(0,0) = vcl_cos(alpha); 
R_y(0,1) = 0;
R_y(0,2) = -vcl_sin(alpha); 
R_y(1,0) = 0;
R_y(1,1) = 1; 
R_y(1,2) = 0; 
R_y(2,0) = vcl_sin(alpha); 
R_y(2,1) = 0; 
R_y(2,2) = vcl_cos(alpha); 

return m_Rotator->Reorient(odf, R_x * R_y); 
}//end OdfRotation


}; //end Class

}//end Functor

template <class TInputImage1, class TInputImage2, class TOutputImage> 
class ITK_EXPORT GeodesicVarianceWithPeakFindingImageFilter :
    public
    BinaryFunctorImageFilter< TInputImage1, TInputImage2, TOutputImage, 
                          Functor::GeodesicVarianceWithPeakFinding< typename TInputImage1::PixelType, 
                                                typename TInputImage2::PixelType, 
                                                typename TOutputImage::PixelType > >
{
public:

  /** Standard class typedefs. */
  typedef GeodesicVarianceWithPeakFindingImageFilter                                             Self;
  typedef BinaryFunctorImageFilter< TInputImage1, TInputImage2, TOutputImage, 
                                    Functor::GeodesicVarianceWithPeakFinding<typename TInputImage1::PixelType, 
                                                         typename TInputImage2::PixelType,
                                                         typename TOutputImage::PixelType >   >                    Superclass;
                           
              
  typedef SmartPointer<Self>                                                  Pointer;
  typedef SmartPointer<const Self>                                            ConstPointer;
  //typedef TOutputImage::RegionType                                            OutputImageRegionType; 

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(GeodesicVarianceWithPeakFindingImageFilter, 
               BinaryFunctorImageFilter);

protected:
  GeodesicVarianceWithPeakFindingImageFilter() {}
  virtual ~GeodesicVarianceWithPeakFindingImageFilter() {}
  
private:

GeodesicVarianceWithPeakFindingImageFilter(const Self&); //purposely not implemented
void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
 

