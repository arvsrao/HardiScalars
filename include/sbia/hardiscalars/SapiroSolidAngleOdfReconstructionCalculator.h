#ifndef __SapiroSolidAngleOdfReconstructionCalculator_h
#define __SapiroSolidAngleOdfReconstructionCalculator_h

#include "itkSymRealSphericalHarmonicRep.h"
#include "string.h"
#include "vnl/vnl_math.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/algo/vnl_svd.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "sphere1000points.h"

/*
* This class implements the constant solid angle ODF reconstruction method of Sapiro, et al.
* For details see their paper: 
*
*    I. Aganj, C. Lenglet, G. Sapiro, E. Yacoub, K. Ugurbil, and N. Harel,  
*    “Reconstruction of the orientation distribution function in single and multiple shell q-ball imaging within constant solid angle,”
*    Magnetic Resonance in Medicine, vol. 64, no. 2, pp. 554–566, 2010. 
*
*/

template< unsigned int SignalSize = 1000, unsigned const short NOrder = 8>

class SolidAngleReconstruction{
public:

typedef itk::SymRealSphericalHarmonicRep<double, NOrder>                        PixelType;
typedef vnl_matrix<double>                                                     MatrixType;
typedef vnl_vector_fixed<double, 1000>                                         VectorType;
typedef vnl_vector_fixed<double, SignalSize >                                  SamplesVecType;
typedef vnl_vector_fixed<double, PixelType::Dimension>                         RshCoefVecType;

SolidAngleReconstruction()
{
  b[0]=0; b[1]=1; b[2] = 6; b[3]=15; b[4]=28;
  A.set_size(PixelType::Dimension, 1000);
  Bmat.set_size(1000, PixelType::Dimension);
  BTmat.set_size(PixelType::Dimension, 1000);
  P.set_size(PixelType::Dimension, PixelType::Dimension);
  L.set_size(PixelType::Dimension, PixelType::Dimension);
  
  for(int i=0; i < 1000; ++i)
    {
         theta.put( i, vcl_acos( points1000Z[i] ) );

         if( points1000Z[i] == 1 || points1000Z[i] == -1 ) //if you're at a pole
         {
            phi.put(i, 0);
         }
         else 
          {
            if( (points1000X[i] == 0) & (points1000Y[i] < 0) )
  	        phi.put(i, double(-1) * vnl_math::pi / double(2) );
            else if( (points1000X[i]==0) & (points1000Y[i] > 0) )
                phi.put(i, vnl_math::pi / double(2));
            else if ( points1000X[i] < 0 )
                phi.put(i, vnl_math::pi + vcl_atan( points1000Y[i] / points1000X[i]) ); 
            else
                phi.put(i, vcl_atan( points1000Y[i] / points1000X[i]) );
         }      
    }

   //create Bmat, BTmat, P, and L matrices of sampled RSH basis functions.
     for(int i=0; i < 1000; ++i)
     {
        
        for(int j = 1; j <= PixelType::Dimension ; ++j)
        { 
            Bmat.put(i, j-1, PixelType::Y(j, theta.get(i), phi.get(i) ) );
            BTmat.put(j-1, i, PixelType::Y(j, theta.get(i), phi.get(i) ) );
        }       
     }
    
     for(int l=0; l <= NOrder; l+=2)
     {   
         //std::cout<< "legendreP(l,0) of l = " << l << " is " << legendreP(l,0) << "\n"; 
        for(int m = -1*l; m <=l ; ++m)
         {
           
           
            std::cout<< l + m + b[ int(double(l)/double(2)) ]<< "\n";
            
            L.put( l + m + b[ int(double(l)/double(2)) ] , l + m + b[ int(double(l)/double(2)) ], double(-1*l*(l+1)) );
            P.put( l + m + b[ int(double(l)/double(2)) ], l + m + b[ int(double(l)/double(2)) ], legendreP(l,0) / double(8 * vnl_math::pi) );   
         }
     
     }
     
     /**
     * Least Squares model fitting is done here. 
     *
     */
     vnl_svd<double> solver(BTmat*Bmat);
     A = P*L* solver.inverse() * BTmat;
     
}//end constructor

virtual ~SolidAngleReconstruction() {}

inline RshCoefVecType ComputeRshCoef(const SamplesVecType samples ) const 
{
      RshCoefVecType shCoef = A*samples;
      shCoef[0]= 1 / (2 * vcl_sqrt(vnl_math::pi));
      return shCoef;
}

private:
  MatrixType                                   A, BTmat, Bmat, P, L; 
  
  VectorType                                   theta, phi, samples; 
  vnl_matrix_fixed<double, 1000, 3>            points;  
  int b[5];  
  
double factorial( int x )
{ 
     double result = 1;
     
     if(x == 0 || x == 1)
        return result;
     else 
     {
        for( int i = int(x) ; i > 0; --i)
        {
           result *= i;
        }
        return result; 
     }
}

double binom( int n, int k)
 {
   double result=1;   
   
   if( n < k)
    {
      std::cerr << "gave invalid pair of integers for combination" << "\n";
      return 0;
    }
   else if( n==k || k==0)
      return 1;
   else
    {   
   
      for( int i = n ; i > n-k ; --i)
        {
           result *= i;
        }
      return ( result / factorial(k) );
    }
 }
 
/* Closed formula of Legendre Function
 * 1/( 2^l * l! ) * sum^{l/2}_{k=0} binom(l,k) * (-1)^k  (2l-2k)! / (l-2k)!  * x^{l-2k} 
 *
 */
double legendreP(int l , double x )
{
     double P_l = 1 /( factorial(l) * vcl_pow( double(2), double(l) ) );  
     double sum = 0;
     
     //Determine (-1)^(l/2)
     int sign = 1;
     if ( vcl_abs(l/2) % 2 == 1) sign = -1;
     
     if (l==0)
         return 1;
     else if( x==0 )
         return sign * ( binom(l, int(double(l)/double(2)))  / vcl_pow(double(2), double(l)) );
     else
      {
        for(int k=0; k <= l/2; ++k)
        {
           sum += sign * binom(l, k)* ( factorial(2*l - 2*k) / factorial(l-2*k)  ) * vcl_pow(double(x), double(l-2*k)); 
        }  
      
       return P_l * sum ;
      }
}

};//end class


#endif
