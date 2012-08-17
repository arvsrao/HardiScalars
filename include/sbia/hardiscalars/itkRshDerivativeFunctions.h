#ifndef __itkRshDerivativeFunctions_h
#define __itkRshDerivativeFunctions_h

#include "itkSymRealSphericalHarmonicRep.h"
#include "itkReplaceSpecialFunctions.h"
#include "spherePoints.h"
#include "itkNumericTraits.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "itkMatrix.h"
#include "itkVector.h"
#include "vector.h"
#include "math.h"

/*
*    itkRshDerivativeunction.h 
*    
*    This file contains functions necessary for calculating first and second derivatives of RSH functions in 
*    the azimuth and elevation directions
*
*
*
*/

namespace itk
{
template< typename TComponentType, unsigned int NOrder = 8>
class ITK_EXPORT RshDerivativeFunctions
{
public:

/** Standard class typedefs. */
typedef RshDerivativeFunctions   Self;
//typedef SmartPointer<Self>                                                           Pointer;
//typedef SmartPointer<const Self>                                                     ConstPointer;

/** Method for creation through the object factory. */
//itkNewMacro(Self);

/** Dimension of the vector space and delcaration of other constants. */

itkStaticConstMacro(numSphPts, unsigned int, 500);
itkStaticConstMacro(numCoeffs, unsigned int, (NOrder+1)*(NOrder+2)/2);

/** Define the component and pixel types */
typedef TComponentType                                                                         ComponentType;
typedef ComponentType                                                                                   OutputPixelType;
typedef itk::SymRealSphericalHarmonicRep< TComponentType, NOrder >                                      RshType;

/** Define vector and matrix types */
typedef itk::Vector< ComponentType, 3 >                                                                 DirectionType;
//typedef itk::Vector< ComponentType, 2 >                                                               VectorType;
typedef vnl_vector< ComponentType >                                                            VectorType; 
typedef vnl_matrix< ComponentType >                                                            MatrixType;

/* Linear Systems Solver  */
typedef itk::SymmetricEigenAnalysis< MatrixType, itk::Vector<ComponentType, 2>, MatrixType>             EigenAnalysisType;

/* angle struct declaration */
struct angles
{
  ComponentType theta; 
  ComponentType phi;
};


/** Default constructor*/
RshDerivativeFunctions()
{

/* Initialize matrix and vector times */

 del_theta_mat.set_size(numSphPts, numCoeffs);
 del_phi_mat.set_size(numSphPts, numCoeffs);
 del_theta2_mat.set_size(numSphPts, numCoeffs);
 del_theta_phi_mat.set_size(numSphPts, numCoeffs);
 del_phi2_mat.set_size(numSphPts, numCoeffs);

  for ( int l = 0; l <= int(NOrder); l+=2)
   {
     for( int m = -l; m <= l; m++)
     {
         for( int i=0; i < numSphPts; ++i ) 
          { 
             dir.SetElement(0, spherePoints[i][0]); 
             dir.SetElement(1, spherePoints[i][1]);
             dir.SetElement(2, spherePoints[i][2]);
             tempAngles  = ConvertSphereCoords( dir );
  
             del_theta2_mat.put(i, getIndex(l,m), del_theta2_Y(l,m, tempAngles.theta, tempAngles.phi));
             del_theta_phi_mat.put(i, getIndex(l,m), del_theta_phi_Y(l,m, tempAngles.theta, tempAngles.phi));
             del_phi2_mat.put(i, getIndex(l,m), del_phi2_Y(l,m, tempAngles.theta, tempAngles.phi));
          }
     }
  
   }
}//end constructor

void GetHessian(RshType rshFunc, MatrixType hessianMat[], int size )
{
  //MatrixType       hessianMat[numSphPts];
  VectorType rshVec(numCoeffs);
    
  for( int i=0; i < numCoeffs; i++) rshVec.put(i, rshFunc[i]);
  
  MatrixType tempMat(2,2);
  VectorType theta2_mat( del_theta2_mat * rshVec);
  VectorType theta_phi_mat( del_theta_phi_mat * rshVec);
  VectorType phi2_mat( del_phi2_mat * rshVec);
  
  
  for(int j=0; j < size; j++)
  {   
      tempMat.put(0,0, theta2_mat.get(j));
      tempMat.put(0,1, theta_phi_mat.get(j));
      tempMat.put(1,0, theta_phi_mat.get(j));
      tempMat.put(1,1, phi2_mat.get(j));

      hessianMat[j] = tempMat; 
     //std::cout<< tempMat[0][0] << " " << tempMat[1][0] << " " << tempMat[0][1] << " " << tempMat[1][1] << "\n";
  }
}

virtual ~RshDerivativeFunctions() {}

static double del_x_LegendreP(int l, int m, double x )
{
  //because l is always even, sign will always equal 0
  int sign_l = 0;
  if ( l%2 == 0 )
     sign_l = 1;
  else 
     sign_l= -1;
  
  
  if ( l == 0 )
      return 0;
  else if ( vcl_abs(m) < l )
      return ((double(m)*x)/(x*x-1)) * LegendreP(l,m,x) - ( double(1) / vcl_sqrt(1-x*x) ) * LegendreP(l,m+1,x);
  else if ( m == -l )
      return ( double(1) / double(doublefactorial(2*l)) ) * double(-1*l) * x * vcl_pow( 1-x*x, (l/2) -1 );
  else 
      return double(sign_l*(-1)*l* doublefactorial(2*l-1)) * x * vcl_pow( 1-x*x, (l/2) -1 );
}

static double del_x2_LegendreP(int l, int m, double x )
{

  int sign_l = 0;
  if ( l%2 == 0 )
     sign_l = 1;
  else 
     sign_l = -1;
  
  if ( l == 0 )
      return 0;
  else if ( vcl_abs(m) < l )
      return (double(1)/ vcl_sqrt(1-x*x)) *( (x/(double(1)-x*x))*LegendreP(l,m+1,x) - del_x_LegendreP(l,m+1,x) ) + (1/(x*x-double(1)))*(m*LegendreP(l,m,x) + double(m-2)*x*del_x_LegendreP(l,m,x) );
  else if ( m == -l )
      return ( double(-1) / double(doublefactorial(2*l)) ) * l * vcl_pow(1-x*x, (l/2)-2) * ( double(1-l)*x*x + double(1));
  else 
      return double( sign_l * -l * doublefactorial(2*l-1) ) * vcl_pow( 1-x*x, (l/2) -2 ) * ( double(1-l)*x*x + double(1));
}

static double del_phi_Y( int l, int m, double theta, double phi )
{
 
  if( m == 0 ) /// Y_l^0
    return 0;
  else if( m < 0 ) /// sqrt2 re(y_l^m)
    return vnl_math::sqrt2 * -m * RshType::K(l,m) * vcl_sin(m*phi) * LegendreP(l,m,vcl_cos(theta));
  else ///(m > 0) sqrt2 im(y_l^m)
    return vnl_math::sqrt2 * RshType::K(l,m) * m * vcl_cos(m*phi) * LegendreP(l,m,vcl_cos(theta));
}

static double del_theta_Y(int l, int m, double theta, double phi )
{

  if( m == 0 ) /// Y_l^0
     return RshType::K(l,0) * -1 * vcl_sin(theta) * del_x_LegendreP(l,m,vcl_cos(theta));
  else if( m < 0 ) /// sqrt2 re(y_l^m)
     return vnl_math::sqrt2 * RshType::K(l,m) * vcl_cos(m*phi) * -1 * vcl_sin(theta) * del_x_LegendreP(l,m,vcl_cos(theta));
  else ///(m > 0) sqrt2 im(y_l^m)
     return vnl_math::sqrt2* RshType::K(l,m) * vcl_sin(m*phi) * -1 * vcl_sin(theta) * del_x_LegendreP(l,m,vcl_cos(theta));
}

static double del_phi2_Y( int l, int m, double theta, double phi )
{
 
  if( m == 0 ) /// Y_l^0
    return 0;
  else if( m < 0 ) /// sqrt2 re(y_l^m)
    return vnl_math::sqrt2 * -m * m * RshType::K(l,m) * vcl_cos(m*phi) * LegendreP(l,m,vcl_cos(theta));
  else ///(m > 0) sqrt2 im(y_l^m)
    return vnl_math::sqrt2* RshType::K(l,m) * -m * m * vcl_sin(m*phi) * LegendreP(l,m,vcl_cos(theta));
}

static double del_theta_phi_Y( int l, int m, double theta, double phi )
{
 
  if( m == 0 ) /// Y_l^0
    return 0;
  else if( m < 0 ) /// sqrt2 re(y_l^m)
    return vnl_math::sqrt2 * -m * RshType::K(l,m) * vcl_sin(m*phi) * -1 * vcl_sin(theta) * del_x_LegendreP(l,m,vcl_cos(theta));
  else ///(m > 0) sqrt2 im(y_l^m)
    return vnl_math::sqrt2* RshType::K(l,m) * m * vcl_cos(m*phi) * -1 * vcl_sin(theta) * del_x_LegendreP(l,m,vcl_cos(theta));
}

static double del_theta2_Y(int l, int m, double theta, double phi )
{
  if( m == 0 ) /// Y_l^0
     return RshType::K(l,0) * -1 * ( vcl_cos(theta)*del_x_LegendreP(l,m,vcl_cos(theta)) - vcl_sin(theta) * vcl_sin(theta) * del_x2_LegendreP(l,m,vcl_cos(theta)) );
  else if( m < 0 ) /// sqrt2 re(y_l^m)
     return vnl_math::sqrt2 * RshType::K(l,m) * -1 * vcl_cos(m*phi) * ( vcl_cos(theta)*del_x_LegendreP(l,m,vcl_cos(theta)) - vcl_sin(theta) * vcl_sin(theta) * del_x2_LegendreP(l,m,vcl_cos(theta)) );
  else ///(m > 0) sqrt2 im(y_l^m)
     return vnl_math::sqrt2* RshType::K(l,m) * -1 * vcl_sin(m*phi) * ( vcl_cos(theta)*del_x_LegendreP(l,m,vcl_cos(theta)) - vcl_sin(theta) * vcl_sin(theta) * del_x2_LegendreP(l,m,vcl_cos(theta)));
}

private:

//private member variables


DirectionType    dir; 
angles           tempAngles;

MatrixType del_theta_mat;
MatrixType del_phi_mat;
MatrixType del_theta2_mat;
MatrixType del_theta_phi_mat;
MatrixType del_phi2_mat;

//private member functions
static unsigned int getIndex(int l, int m)
{
//int b[0]={0,1,6,15,28}; 
  switch(l)
  {
     case 0:
        return l+m+0;
        break;
     case 2:
        return l+m+1;
        break;
     case 4:
        return l+m+6;
        break;
     case 6:
        return l+m+15;
        break; 
     case 8:
        return l+m+28;
        break;
     default:
        std::cerr<< "ERROR: l > 8 or odd" << std::endl;
        return 0;
  }
}

static unsigned int factorial(unsigned int x)
{
    unsigned int value = 1;

    for(unsigned int i = 2; i <= x; i++)
    {
        value = value * i;
    }

    return value;
}

static int doublefactorial(unsigned int n)
{
    int i;
    int res = 1;      
    for(i=n;i>=1;i-=2)      
     {
       res *=i;
     }
    return res;
}

static angles ConvertSphereCoords( DirectionType dir )
{
   angles temp; 
   if( ( dir.GetElement(2) == 1) || ( dir.GetElement(2) == -1) ) //if you're at a pole
   {
       temp.theta = 0;
       temp.phi    = vcl_atan2(dir.GetElement(1) , dir.GetElement(0));
   }
   else 
   {
       temp.theta  = vcl_acos( dir.GetElement(2)); 
       temp.phi    = vcl_atan2(dir.GetElement(1) , dir.GetElement(0));
   }
     
return temp;
} 

};//end class

}//end namespace
#endif
