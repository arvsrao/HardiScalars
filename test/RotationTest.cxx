#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNiftiImageIO.h>
#include <itkSymRealSphericalHarmonicRep.h>
#include "string.h"
#include "math.h"

#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include <itkRSHPixelReorientationOperator.h>
#include "itkPeakFindingCalculator.h"
#include <itkMatrix.h>
#include <itkVector.h>
#include "sphere1000points.h"
/****************************************************
                       Typedef's and Constants
 ****************************************************/

unsigned const short Dimension = 3;
unsigned const short NOrder = 8; 

//Pixel Types
typedef itk::SymRealSphericalHarmonicRep<double, NOrder> PixelType;

//Image Types
typedef itk::Image<PixelType, Dimension>             ImageType;

//IO types
typedef itk::ImageFileReader<ImageType>              ReaderType; 


/****************************************************
                       MAIN
 ****************************************************/


template<class TInput>
class Rotation
{
public:
  
  typedef TInput                                                   PixelType; 
  
  /** Container to hold peak directions of the each ODF */
  typedef vnl_vector_fixed<double, 3>                              PeakDirectionType;  
  typedef itk::PixelReorientationOperator<PixelType>               RotationOperatorType;
  typedef typename   RotationOperatorType::Pointer	           RotationOperatorPointer;
  
  
  typedef itk::PeakFindingCalculator< PixelType, double, double>   PeakFindingCalculator;
  typedef typename   PeakFindingCalculator::Pointer                PeakFindingCalculatorPointer;
  typedef typename   PeakFindingCalculator::InputType              PeakFindingInputType;                 ;
  
  typedef itk::Matrix<double, 3, 3>                                MatrixType; 
  typedef vnl_matrix<double>                                       VnlMatrixType;
  typedef vnl_vector_fixed<double, 10000>                          VectorType;
  typedef vnl_vector<double>                                       SamplesVecType;
  int                                                              count;
    
  
Rotation() 
{
    count =0;
    m_Rotator = RotationOperatorType::New();
    m_pFinder = PeakFindingCalculator::New();
    mat.set_size(10000, TInput::Dimension);

    for(int i=0; i < 10000; ++i)
    {
       points.put(i,0, points10000X[i]);
       points.put(i,1, points10000Y[i]);
       points.put(i,2, points10000Z[i]);
    }

   for(int i=0; i < 10000; ++i)
    {
         theta.put( i, vcl_acos( points10000Z[i] ) );

         if( points10000Z[i] == 1 || points10000Z[i] == -1 ) //if you're at a pole
         {
            phi.put(i, 0);
         }
         else 
          {
            phi.put(i, atan2( points10000Y[i], points10000X[i] ));
         }      
    }

   //create mat matrix of sampled SH basis functions and fill sine and cosine vectors;
     for(int i=0; i < 10000; ++i)
     {
        sine.put(i, vcl_sin( theta.get(i) ) );
        cosine.put(i, vcl_cos( theta.get(i) ) );
        
        for(int j = 1; j <= TInput::Dimension ; ++j)
        { 
            mat.put(i, j-1, TInput::Y(j, theta.get(i), phi.get(i) ) );
        }       
     }
     std::cerr<< " finished with constructor" << "\n";
}

virtual ~Rotation() {}


inline PixelType RotationFunction( const PixelType A ) 
{
    
  vnl_vector_fixed<double , TInput::Dimension>      input;
  // std::cerr << "Input to GetVariance Operators" << std::endl;
  // std::cerr << A << std::endl << std::endl;
   // instantiate input vnl_vector_fixed of
      for(int i=0; i < TInput::Dimension ; ++i)
       {
          input.put( i , A[i] );
       }
   
  
   if( input.one_norm() == 0 )
       return A;
   else
   {

/*
 * Peakfinding and ODF rotation is done here
 *
 */     
     OneNorm = 0;
     sum=0;
     samples  = mat * input;
     int                p;
     PixelType          C;
     

     for(int i=0; i <10000; ++i)
      {
         if ( sum < samples.get(i) )
	  { 
	     sum = samples.get(i);
	     p=i;
          }
      }

     peak.put(0, points10000X[p]);
     peak.put(1, points10000Y[p]);
     peak.put(2, points10000Z[p]);
         
     
     std::cout<< "The peak (locating the max point) is: ("<< peak.get(0) << ", " << peak.get(1) << ", "<< peak.get(2) << ") "<< "with max:  "<< sum <<"\n\n";
    
     std::vector< double >  odPeaks; 
     PeakFindingInputType B;
     
     for(int j=0; j < 45; ++j) B[j] = A[j];
     
     m_pFinder->SetCoefficients( B );
     m_pFinder->SetPeakFiltering( 1 );
     m_pFinder->SetAcceleration( 0 );
    
     int nmPeaks = m_pFinder->GetRawPeaks( odPeaks );
     
     std::cout<< "(Using peak finding code) peak is: ("<< odPeaks[0] << ", " << odPeaks[1] << ", "<< odPeaks[2] << ")  " << " with max: " << odPeaks[3] <<"\n\n";
   
     peak.put(0, odPeaks[0]);
     peak.put(1, odPeaks[1]);
     peak.put(2, odPeaks[2]);
      
    //Print out peak of rotated ODF
     C = OdfRotation(A, peak);
    
     for(int j=0; j < 45; ++j) B[j] = C[j];
     odPeaks.clear();
     m_pFinder->SetCoefficients( B );
     m_pFinder->SetPeakFiltering( 1 );
     m_pFinder->SetAcceleration( 0 );
    
     nmPeaks = m_pFinder->GetRawPeaks( odPeaks );
     
     std::cout<< "(Using peak finding code) peak of rotated Odf is: ("<< odPeaks[0] << ", " << odPeaks[1] << ", "<< odPeaks[2] << ")\n\n";

     if( vcl_acos( peak.get(2) ) <.0001 )
        return A;
     else 
        return C;
    }
}

/*Private Member functions*/

  


//Creates a rotation matrix that takes p to q.
inline PixelType OdfRotation(PixelType odf, PeakDirectionType peak_temp)
{

//Compute Euler Axis;
float alpha = vcl_acos(peak_temp.get(2));
float norm_axis = vcl_sin(alpha);
PeakDirectionType euler_axis(-1*peak_temp.get(1)/norm_axis, peak_temp.get(0)/norm_axis, 0);

//Compute rotation matrix
MatrixType A, RotMat;

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

std::cout << A * RotMat *A << "\n\n";

std::cout << A * RotMat * A * peak_temp << "\n\n";

std::cout << "angle of peak with North pole is: " << 180 * alpha / vnl_math::pi <<"\n\n";


return m_Rotator->Reorient(odf, A*RotMat*A); 
}//end OdfRotation

private:

PeakDirectionType                               peak;
double                                          OneNorm, sum;
VnlMatrixType                                   mat; 
VectorType                                      theta, phi, samples, sine, cosine; 
vnl_matrix_fixed<double, 10000, 3>              points;

RotationOperatorPointer m_Rotator;
PeakFindingCalculatorPointer m_pFinder;

}; //end Class



int main(int argc, char* argv[])
{


if( argc < 2 )
  return EXIT_FAILURE;

const int k = atoi(argv[1]); 
    

double CoeffMat[2][45] = 
{
// 56 50 44
{0.28209481, 0.044422667, -0.00029324778, -0.079764642, -0.018271139, 0.17813154, -0.06166628, -0.012195878, -0.011008681, 0.013260418, 0.031245476, 0.012355535, -0.048984621, -0.0078436024, 0.035731435, -0.0086415056, -0.00090472057, 0.012633592, 0.0065435679, 0.0032198436, -0.0071723103, -0.00326701, -0.0084158853, 0.003511135, 0.0075718607, -0.008445885, 0.00069137174, -0.013899945, 0.0023685631, 0.0013487741, 0.0069391346, 0.0006436322, 0.0072214734, -0.0061521847, 0.0029057232, -0.0029879399, -0.00027989462, -0.0010987743, -0.003424343, 0.0011984621, 0.0011805261, 0.0001017326, 0.0021068042, -0.0040239459, -0.001655195},

// 73 50 44 
{0.28209481, 0.031946689, 0.012138316, -0.037802678, -0.0062907366, -0.064339511, -0.0092263399, 0.0076692216, -0.0070908326, -0.0015781252, 0.012303672, 0.010785273, 0.014094436, 0.0020881153, -0.018144067, -0.0013923909, 0.0018931613, 0.003348592, -0.00076843746, 0.0011416436, -0.00021286854, -7.1700553e-05, 0.00024459604, -0.00021596726, -0.002409227, -0.00059268484, 0.0019698753, 0.00068739493, -0.0003986902, 0.00097694353, -0.0010657698, -0.00075136009, -0.00041485927, 0.0018464801, 0.0027769054, 0.0031351894, 0.00093688653, 0.0013665569, 0.0019333645, 0.00019473411, -0.0011095481, -0.0016485921, 0.00050739496, 0.000482973, 0.0011656274}
};

PixelType RshPixel;
for(int i =0; i <45; ++i) RshPixel[i] = CoeffMat[k-1][i];

Rotation<PixelType> i_Rotator; 
PixelType C = i_Rotator.RotationFunction( RshPixel );

for(int i =0; i <45; ++i) std::cout << C[i] << ", ";
std::cout<< "\n"<<std::endl;

return 0;
}//end main

