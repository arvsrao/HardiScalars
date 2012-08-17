#define VERBOSE_CXX
#define INFO_FILE

#include "itkSymRealSphericalHarmonicRep.h"
#include "itkGeodesicConcentrationFilter.h"
#include "itkRSHPixelReorientationOperator.h"
#include "itkRshDerivativeFunctions.h" 
#include "itkSymmetricEigenAnalysis.h"
#include "itkVector.h"
#include "string.h"
//#include "tclap/CmdLine.h"
#include "iostream.h"
#include "math.h"
#include "spherePoints.h"

#include "itkNumericTraits.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_matrix.h"

/*
 * ComputeGCatSpherePoints.cxx
 *  
 * This executable computes Guassian Curvature of RSH functions at 500 equally spaced
 * sphere points defined in "spherePoints.h". As of 6/25/2012 it is in testing. 
 *
 */

//Template parameters are only handled at compile time. 
float shRepCoeffs[4][45] =
          {
            {
               2.82094806e-01,  4.59387809e-01,  7.62639241e-03, -2.65089244e-01,
              -1.70459377e-03, -2.27291649e-03,  3.38489264e-01,  8.24082829e-03,
              -2.55495191e-01, -8.35784525e-03,  1.70923799e-01,  1.99015043e-03,
               1.02482527e-03, -5.59562538e-03, -3.80363222e-03,  1.70554623e-01,
               5.64062875e-03, -1.26915440e-01, -5.26009873e-03,  1.14953570e-01,
               4.45495034e-03, -7.90442452e-02, -1.32022612e-03, -3.45894368e-04,
               3.95073835e-03,  9.47112043e-04, -6.37689978e-03, -3.68007110e-03,
               4.89506535e-02,  2.26373016e-03, -3.73514183e-02, -2.18456332e-03,
               3.32550853e-02,  1.17909722e-03, -3.10148522e-02, -6.25214074e-04,
               2.12846715e-02,  3.60342383e-04,  3.37901074e-05, -1.33388466e-03,
              -1.39679192e-04,  2.47698487e-03,  2.21211842e-04, -3.72869126e-03,
              -2.16668053e-03
            },

            {
               2.82094806e-01, -4.56637770e-01, -2.07055407e-03, -2.66447693e-01,
              -3.83342733e-03, -2.03245622e-03,  3.29349637e-01,  6.67131925e-03,
               2.56829888e-01,  2.59662815e-03,  1.74307019e-01,  4.18567704e-03,
               2.51308928e-04,  3.66590638e-03,  3.40090320e-03, -1.61080286e-01,
              -6.84451079e-03, -1.24465309e-01, -5.27192699e-03, -1.18184127e-01,
              -1.89521688e-03, -8.28848481e-02, -2.20265985e-03,  4.79890703e-04,
              -2.14970903e-03,  4.34971036e-04, -1.83163770e-03, -3.47044342e-03,
               4.55690846e-02,  3.02502280e-03,  3.43799517e-02,  3.41011304e-03,
               3.30946259e-02,  2.10025162e-03,  3.35046090e-02,  6.53992000e-04,
               2.37989724e-02,  3.83429928e-04, -1.72163273e-04,  3.30030743e-04,
              -7.50337844e-04,  3.65657092e-04, -9.34689073e-04,  3.83923936e-04,
               2.20232457e-03
            },

            {
               2.82094806e-01, -1.84499600e-03, -4.37358674e-03, -2.57216394e-01,
              -2.21597566e-03, -5.69619704e-03,  3.87160391e-01, -5.81498304e-03,
              -4.88071050e-03,  3.55852279e-03,  1.55444667e-01,  1.25033641e-03,
               2.96612363e-03,  2.67166062e-03, -6.83821831e-03, -1.60360534e-03,
              -3.04495008e-03, -1.28530160e-01,  2.44779466e-03,  5.92144812e-03,
              -4.25606035e-04, -6.48434013e-02,  5.08837402e-04, -1.10527303e-03,
              -2.69610246e-05,  2.88396259e-03, -2.25038873e-03, -7.32504530e-03,
               1.14771903e-01, -3.58927739e-03, -2.93563283e-03,  1.99614253e-04,
               2.79063080e-02,  2.29189231e-04, -2.89570610e-03, -1.19774416e-03,
               1.54047897e-02, -7.50655425e-04,  2.04480690e-04, -1.01385545e-03,
              -9.12652526e-04, -5.92372544e-06,  1.25703122e-03,  2.75207148e-03,
              -3.13866278e-03
            },
            
            {
              .222315139564, -.166693234407, 0.851766149229e-2, .405768890937, 
              .234067800088, 0.43792184354e-3, 0.29281169775e-1, 0.224117918194e-3, 
              -.112278771026, 0.509434374726e-2, .12893490089, 0.433954824817e-1, 
              -0.779228726611e-3, -0.829873209978e-1, -0.943169153778e-4, 0.108872811721e-1, 
              0.83592673924e-2, 0.950815589252e-1, 0.224815789131e-2, -0.486475915346e-1, 
              0.191965115135e-1, .183717027584, 0.624021906204e-2, -0.204536028673e-2, 
             -.166080403816, 0.368281318914e-2, 0.416128659131e-1, 0.148847770717e-1, 
             0.216923815881e-1, -0.116101172405e-1, -0.60335880558e-1, 0.485929634959e-3, 
             0.288728364256e-1, -0.683209490229e-3, 0.915442739555e-1, 0.98380006569e-2, 
             .208472701934, -0.547559498446e-1, 0.548132769607e-2, 0.305184665611e-1, 
             -0.744708012157e-2, 0.8726741661e-1, -0.192179113597e-2, -0.10700691283e-1, 
             0.518404877475e-2
             
            }
          };

const unsigned int    NOrder = 8;
const unsigned int    Dimension = 3;
const int             arrayLength = 500;

typedef float         ComponentType;

typedef itk::SymRealSphericalHarmonicRep< ComponentType, NOrder >                                       RshPixelType;
typedef itk::RshDerivativeFunctions< ComponentType, NOrder>                                              FunctionsClass; 
typedef ComponentType                                                                                   OutputPixelType;
typedef itk::Vector< ComponentType, Dimension >                                                         DirectionType;
//typedef itk::Functor::GeodesicConcentration< RshPixelType, DirectionType, OutputPixelType >             CalcType;
typedef FunctionsClass::MatrixType                                                                      MatrixType;
typedef itk::SymmetricEigenAnalysis< MatrixType, itk::Vector<ComponentType, 2>, MatrixType>             EigenAnalysisType;

typedef itk::PixelReorientationOperator< RshPixelType >                                                 RotationOperatorType;
typedef RotationOperatorType::Pointer	                                                                RotationOperatorPointer;

struct angles
{
    ComponentType theta;
    ComponentType phi;
};

int main(int argc, char* argv[])
{

FunctionsClass derivatives;
RshPixelType   shCoefs[ 4 ];
MatrixType what[2];
float tempa[4] = {1,1,1,1};


itk::Vector<ComponentType, 2>         vec;
EigenAnalysisType eigenAnalysis(2);

//declare temp vnl_matrix;
MatrixType hessArray[arrayLength];
for ( int i = 0; i < 4; i++ ) shCoefs[i] = shRepCoeffs[i];

derivatives.GetHessian(shCoefs[3], hessArray, arrayLength);

for( int j=0; j < arrayLength; ++j)
{
   eigenAnalysis.ComputeEigenValues( hessArray[j], vec);
   if( (vec[0] <0 & vec[1] > 0) | (vec[0] > 0 & vec[1] < 0) )
       std::cout<<" saddle point @ " << j << "\n";
}
    
return EXIT_SUCCESS;
}//end main

//for( int i=0; i < arrayLength; ++i ) 
//{ 
//  dir.SetElement(0, spherePoints[i][0]); 
//  dir.SetElement(1, spherePoints[i][1]);
//  dir.SetElement(2, spherePoints[i][2]);
//  tempAngles  = ConvertSphereCoords( dir );
//  
//  mat.put(0,0, derivatives.del_theta2_Y(4,3, tempAngles.theta, tempAngles.phi));
//  mat.put(0,1, derivatives.del_theta_phi_Y(4,3, tempAngles.theta, tempAngles.phi));
//  mat.put(1,0, derivatives.del_theta_phi_Y(4,3, tempAngles.theta, tempAngles.phi));
//  mat.put(1,1, derivatives.del_phi2_Y(4,3, tempAngles.theta, tempAngles.phi));
//  eigenAnalysis.ComputeEigenValues(mat, vec);
//  std::cout<<  vec  << "\n"; 
//}


//    #ifdef INFO_FILE
//        FILE *fd = fopen("/sbia/home/raoarvin/GcValues.txt", "w");
//        //fprintf( fd, "\n\n\nRotated Coefficients Are:   \n\n");    
//        for(int i=0; i< arrayLength; ++i)
//        {  
//           //std::cout << outputGC[i] << "\n";
//           fprintf( fd, "%.12g \n", outputGC[i]);
//        }
// 
//        fclose(fd);
//    #endif

