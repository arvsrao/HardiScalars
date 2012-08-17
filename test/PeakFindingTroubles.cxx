#include "itkSymRealSphericalHarmonicRep.h"
#include "itkPeakFindingCalculator.h"
#include "itkVector.h"
#include "itkMatrix.h"
#include "itkRSHPixelReorientationOperator.h"
#include "vnl/vnl_matrix.h"
#include "stdio.h"

#define INFO_FILE

const int    shrMaxOrder = 8;
const int    numberVoxls = 2;

const int    numShCoeffs = ( shrMaxOrder + 1 ) * ( shrMaxOrder + 2 ) / 2;

typedef itk::SymRealSphericalHarmonicRep< double , 8 >     PixelType;
typedef itk::PeakFindingCalculator< PixelType, double, double >            PeakFindingCalculator;

typedef itk::PixelReorientationOperator< PixelType >                       RotationOperatorType;

typedef itk::Matrix<double, 3, 3>                                          MatrixType; 
typedef vnl_matrix<double>                                                 VnlMatrixType;

int testPeakFinding()
{
  
    double  shRepCoeffs[numberVoxls][numShCoeffs] =
    { 
        // PixelOne
        {.282094806, 0.586144924e-1, -0.606478080e-1, .107649803, -0.403662883e-1, -0.617511384e-1, -0.170150511e-1, 0.118303485e-2, 0.182855856e-1, -0.303145382e-2, 0.280970745e-1, -0.650472520e-2, -0.247074943e-1, -0.444193780e-1, -0.326748565e-1, -0.790867954e-2, 0.534386933e-2, -0.702696573e-2, 0.930769966e-4, -0.111849385e-2, 0.699318759e-2, 0.338856247e-2, 0.304164458e-2, -0.400963752e-2, 0.859028252e-3, -0.809718482e-2, -0.287311687e-2, -0.801092573e-2, -0.395215768e-2, 0.377272605e-2, 0.616312609e-3, -0.110100827e-2, -0.662258390e-4, 0.872983132e-3, -0.820796646e-3, 0.544747629e-3, 0.128666428e-2, 0.130909891e-2, -0.272324891e-3, 0.249228696e-2, -0.124018145e-3, 0.220012176e-2, -0.178599890e-3, 0.228770450e-2, 0.266295206e-2},
        
        //PixelTwo
        {.282094806, 0.188417975e-1, 0.559826121e-1, .133148223, 0.194038078e-1, 0.509412736e-1, -0.226358995e-1, 0.562748406e-2, -0.512774495e-4, 0.187420230e-1, 0.330971852e-1, 0.155681847e-1, 0.286710821e-1, -0.336171985e-1, 0.149136288e-1, -0.684208469e-2, -0.186275179e-2, -0.346410507e-2, 0.109557156e-2, 0.142522086e-2, 0.744172605e-2, 0.219848170e-2, 0.698357774e-2, 0.284683588e-2, -0.234303949e-2, 0.425443659e-2, 0.274209515e-2, -0.406520441e-2, 0.225594238e-3, -0.123606075e-2, -0.109349284e-2, -0.170162297e-2, 0.258140429e-2, 0.718621013e-3, -0.211267290e-2, 0.209630397e-2, -0.197931458e-4, 0.224253861e-2, -0.240332773e-2, -0.641495455e-3, 0.239522406e-2, 0.694274495e-3, 0.169952726e-2, -0.283206184e-2, -0.205701659e-2}      
     };
        
    int i,j;
    std::vector< double >            odPeaks; // each peak by 4 components: vx, vy, vz, mag
    odPeaks.clear();

    PeakFindingCalculator::Pointer       pFinder =  PeakFindingCalculator::New(); 
    PeakFindingCalculator::InputType     shCoefs[numberVoxls];
    for(i=0; i<numberVoxls; ++i) shCoefs[i] = shRepCoeffs[i];

           
    // ------------------------------- the test -------------------------------
        
        pFinder->SetCoefficients( shCoefs[1] );
        pFinder->SetPeakFiltering( 1 );
        pFinder->SetAcceleration( 0 );
        
        int   nmPeaks   = pFinder->GetRawPeaks( odPeaks );
                
        for ( j=0; j<36; ++j) 
        { 
            if( j%4 == 3 )
             {
                std::cout<<"\nAmplitude of peak: "<< (j+1)/float(4)<< " is: "<< odPeaks[j] <<"\n";
                std::cout<<"...with peak: ("<< odPeaks[j-3] << ", "<< odPeaks[j-2] << ", "<< odPeaks[j-1] << ") "<<"\n"; 
             }  
        }
   
   #ifdef INFO_FILE
    FILE *fd = fopen("PeakTestingInfo.txt", "w");
    
    for(i=0; i<36; ++i)
     {
        if( i%4 == 3 )
            {
                fprintf(fd, "Amplitude of peak ");
                fprintf(fd, "%.1g", (i+1)/double(4) );
                fprintf(fd, " is ");
                fprintf(fd, "%.12g, ",odPeaks[i]);
                fprintf(fd, "|   (%.12g, ",odPeaks[i-3]);
                fprintf(fd, "%.12g, ",odPeaks[i-2]);
                fprintf(fd, "%.12g )\n",odPeaks[i-1]);
            }  
     }
   #endif 
    
      
    double RotMat[numberVoxls][9] = { {0.328421035298684050, 0.600766210416498803,-0.728848121399999993, 
                                    -0.877446746800000010, 0.479674062900000031, 0, 
                                     0.349609539628970456, 0.639525413033721457, 0.684675409200000029 },
        
                                    { -.133661593285519770, -.974379218673503300, -.180885921899999996,
                                       .990722086999999974,  -.135903445099999992,  0, 
                                     -0.245830199562995372e-1, -.179207678053687008, .983504084000000001}
                                    };
    
    
    VnlMatrixType VnlRotArray[numberVoxls];
    MatrixType ItkRotArray[numberVoxls];
          
    for(i=0; i<numberVoxls; ++i) 
    {
       VnlRotArray[i].copy_in( RotMat[i] );
       ItkRotArray[i] = VnlRotArray[i] ; 
    }
   
    RotationOperatorType::Pointer   RotationOperator =  RotationOperatorType::New();
    
    PixelType odf_out, odf_in;
    for(i=0; i< numShCoeffs; ++i) odf_in[i] = shCoefs[k][i];
    
    odf_out = RotationOperator->Reorient( odf_in, ItkRotArray[k] );
 
 
    #ifdef INFO_FILE
    
        fprintf( fd, "\n\n\nRotated Coefficients Are:   \n\n");    
        for(i=0; i< numShCoeffs; ++i)
        {
           fprintf( fd, "%.12g, " , odf_out[i]);
        }
   
        fclose(fd);
    #endif
  
    return 0;
}

using namespace itk;
int main( int argc, char * argv[] )
{
const int k = atoi( argv[1] );

    if( argc < 2 )
       {
         std::cout<< "Missing index!!"<<"\n";
         return EXIT_FAILURE;
       }
    else
       { 
         testPeakFinding( k );
         return EXIT_SUCCESS;
       }
}
