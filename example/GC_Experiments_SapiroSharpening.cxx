#include "itkSymRealSphericalHarmonicRep.h"
#include "SapiroSolidAngleOdfReconstructionCalculator.h"
#include "string.h"
#include "stdio.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_matrix_fixed.h"
#include "SpherePoints500.h"

unsigned const short NOrder = 8; 
unsigned const short Dimension = (1+NOrder) * (NOrder+2)/2;
typedef SolidAngleReconstruction<500, NOrder>                                  CalculatorType;
typedef vnl_matrix<double>                                                     MatrixType;
typedef vnl_vector_fixed<double, 500>                                          VectorType;
double delta = .001;

double SmoothingFunction( double x )
{
    
    if( x < 0)
        return delta/2;
    else if( x < delta & x >=0)
        return delta/double(2) + x*x/double(2*delta);
    else if( x >=delta & x < 1-delta )
        return x; 
    else if( x >= (1-delta) & x <1 )
        return 1- delta/double(2)  - (1-x)*(1-x)/double(2*delta) ;
    else 
        return 1 - delta/double(2);   
}

int main(int argc, char* argv[])
{
    
    /* Print RSH coefficents to file 
     * Print Parameters coefficients to file
     *
     */
    FILE *fd, *fp; 
    fd = fopen("/Users/arvind/Desktop/Coefficents_SapiroSharpened_Odfs.h", "w");
    fp = fopen("/Users/arvind/Desktop/Parameters_GC_Sapiro.h", "w");
    
    if( fp == NULL )
    {
        std::cerr << "failed to open Parameters file "<<"\n";
        return EXIT_FAILURE;
    }
    
    if( fd == NULL)
    {
        std::cerr << "failed to open Coefficients file "<<"\n";
        return EXIT_FAILURE;
    }
    
    
    
    MatrixType                                points, Tensor_one, Tensor_two, U, W, V, A, B;
    VectorType                                samples;
    double                                    c, d, E_1, E_2, det_tensor_one, det_tensor_two;
    CalculatorType                            calc;  
    vnl_vector_fixed<double, Dimension>       shCoef, rotatedShCoef; 
    vnl_vector_fixed<double, NOrder>          parameters;
    
    points.set_size(500,3);
    Tensor_one.set_size(3,3);
    Tensor_two.set_size(3,3); 
    U.set_size(3,3);
    W.set_size(3,3);
    V.set_size(3,3);
    A.set_size(3,3);
    B.set_size(3,3);
    
    U.fill(0); W.fill(0); V.fill(0); Tensor_one.fill(0); //fill V, U, W with zeros
    U.put(1,1, 1); V.put(1,1, 1); Tensor_one.put(2,2,10); 
    
    double tensor_array[9] = { 1,0,0,
                               0,1,0,
                               0,0,10 };
    Tensor_two.copy_in(tensor_array);
    
    for(int i=0; i < 500; ++i)
    {
        points.put(i,0, points500X[i]);
        points.put(i,1, points500Y[i]);
        points.put(i,2, points500Z[i]);
    }
    
    
    for (double beta=1; beta <=10; ++beta) 
    {
        std::cout<< " beta = " << beta << "\n";
        for(double alpha = 0; alpha < vnl_math::pi/double(2); alpha += vnl_math::pi/double(20) )
        {
            
          // std::cout<< " alpha = " << alpha << "\n";
            for(double a=.1; a<=1; a+=.1)
           {
          // std::cout<< " a = " << a << "\n";
 
                Tensor_one.put(0,0, beta) ;
                Tensor_one.put(1,1, beta) ;
                
                U.put(0,0, vcl_cos(alpha));
                U.put(0,2, -1 * vcl_sin(alpha));
                U.put(2,0, vcl_sin(alpha));
                U.put(2,2, vcl_cos(alpha));
                
                
                V.put(0,0, vcl_cos(vnl_math::pi/double(2) - alpha));
                V.put(0,2, -1 * vcl_sin(vnl_math::pi/double(2) - alpha));
                V.put(2,0, vcl_sin(vnl_math::pi/double(2) - alpha));
                V.put(2,2, vcl_cos(vnl_math::pi/double(2) - alpha));
                
                
                det_tensor_one = vcl_sqrt(vnl_determinant(Tensor_one));
                det_tensor_two = vcl_sqrt(vnl_determinant(Tensor_two));
                
                W = U.transpose() * Tensor_two * U;
                
                c = det_tensor_one /(det_tensor_one + det_tensor_two);
                d = det_tensor_two /(det_tensor_one + det_tensor_two);
                
                //Signal Creation  
                for( int i = 0; i < 500 ;  ++i)
                {
                    E_1 = a*(det_tensor_one / vcl_pow(double(2*vnl_math::pi), double(1.5) ) ) * vcl_exp( double(-.5) * dot_product( points.get_row(i) , Tensor_one * points.get_row(i) ) ) ;
                    E_2 = (1-a)*(det_tensor_two / vcl_pow(double(2*vnl_math::pi), double(1.5) ) ) * vcl_exp( double(-.5) * dot_product( points.get_row(i) , W * points.get_row(i) ) ) ;
                    
                    samples.put(i, vcl_log( -1 * vcl_log( SmoothingFunction( E_1 + E_2 ))) );
                  // samples.put(i, vcl_log( -1 * vcl_log( SmoothingFunction( E_1))) );
                } 
                
                shCoef = calc.ComputeRshCoef(samples);
                
                //Write to File Coefficients shCoef
                fprintf(fd, "{ ");
                for(int i =0; i < Dimension ; ++i)
                { 
                    if(i==Dimension-1)
                        fprintf(fd, "%.15f},\n ", shCoef.get(i));
                    else 
                        fprintf(fd, "%.15f, ", shCoef.get(i)); 
                }
                
                
                A = V.transpose() * Tensor_one * V;
                B = V.transpose() * W * V;
                
                //Rotated Signal Creation  
                for( int i = 0; i < 500;  ++i)
                {
                    E_1 = a * (det_tensor_one / vcl_pow(double(2*vnl_math::pi), double(1.5) ) ) * vcl_exp( double(-.5) * dot_product( points.get_row(i) , A * points.get_row(i) ) ) ;
                    E_2 = (1-a) * (det_tensor_two / vcl_pow(double(2*vnl_math::pi), double(1.5) ) ) * vcl_exp( double(-.5) * dot_product( points.get_row(i) , B * points.get_row(i) ) ) ;
                    
                    samples.put(i, vcl_log( -1 * vcl_log( SmoothingFunction( E_1 + E_2) )) );
                }   
                
                rotatedShCoef = calc.ComputeRshCoef(samples); 
                
                parameters[0] = beta;                                 // width of the lobe being varied
                parameters[1] = 0; //calc.ComputeGC( shCoef );        // GC of [0,0,1] peak    
                parameters[2] = 0; //calc.ComputeGC( rotatedShCoef);  // GC of [1,0,0] peak
                parameters[3] = 0; //calc.ComputeGFA(shCoef);         // GFA of whole ODF
                parameters[4] = 90 * ( vnl_math::pi - 2 * alpha ) / vnl_math::pi ; //angle between peaks in degrees
                parameters[5] = a ;                                   //   partial volume of [0,0,1] peak
                parameters[6] = 1-a;                                  //   partial voluem of [1,0,0] peak
                
                //Write parameters to file
                for(int i =0; i < 7 ; ++i)
                { 
                    if(i==7-1)
                        fprintf(fp, "%.15f\n ", parameters.get(i) );
                    else 
                        fprintf(fp, "%.15f  ", parameters.get(i) ); 
                }
                
                
            }
        }
    }
    
    

    fclose(fd); 
    fclose(fp);
    
}//end main
