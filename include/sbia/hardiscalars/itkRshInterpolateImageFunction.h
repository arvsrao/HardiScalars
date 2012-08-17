#ifndef __itkRshInterpolateImageFunction_h
#define __itkRshInterpolateImageFunction_h

#define VERBOSE 

#include "itkConjugateGradientOptimizer.h"
#include "itkSingleValuedCostFunction.h"
#include "itkImageAdaptor.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkNthElementImageAdaptor.h"
#include "itkMacro.h"

#include "itkVector.h"
#include "itkMatrix.h"
#include "math.h"
#include "spherePoints.h"
#include "vector.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_math.h"
#include "vnl/vnl_matrix.h"

/*
*  
*  This class implements a method of BSpline interpolation for RSH functions
*  represented in the real spherical harmonic basis. The interpolation is
*  applied coefficient-wise. As opposed to orientation-wise. 
*   
*  Coefficients are found using the conjugate gradient method of finding the minimal
*  solution to a Linear system:
*
*
*/


class conjugateCostFunction : public itk::SingleValuedCostFunction 
{
public:
    
    typedef conjugateCostFunction                    Self;
    typedef itk::SingleValuedCostFunction            Superclass;
    typedef itk::SmartPointer<Self>                  Pointer;
    typedef itk::SmartPointer<const Self>            ConstPointer;
    
    itkNewMacro( Self );
    itkTypeMacro( conjugateCostFunction, SingleValuedCostFunction );
    
    typedef Superclass::ParametersType              ParametersType;
    typedef Superclass::DerivativeType              DerivativeType;
    
    typedef vnl_vector<double>                      VectorType;
    typedef vnl_matrix<double>                      MatrixType;
    
    typedef double MeasureType ;
    
    //Default Constructor
    conjugateCostFunction(){}
        
    conjugateCostFunction( MatrixType input , VectorType output )
    {
        A = input;
        F = output;
       
        m_Dimension = output.size(); 
        
        //create column vector from samples of F.
        F_mat.set_size( m_Dimension, 1);
        for( int i = 0; i < m_Dimension; ++i) F_mat.put(i,1, output[i] );       
    }
    
    /** 
     *    For solving linear systems of the type
     *
     *           A*x = b
     * 
     *  The objective function is the quadratic form:
     *
     *  1/2 x^T A x - x^T b
     *
     *  Where A is represented as VNL MATRIX and 
     *  b is represented as a VNL VECTOR
     *
     *  An excample system:
     *
     *     | 3  2 ||x|    | 2| 
     *     | 2  6 ||y| =  |-8| 
     *
     *
     *  the solution is the vector | 2 -2 |
     *
     */ 
    double GetValue( const ParametersType & position ) const
    { 
       MatrixType vec(m_Dimension,1);
       for( int i = 0; i < m_Dimension; ++i) vec.put(i,1, position[i] );
       
       MatrixType val = vec.transpose() * ( .5 * A * vec -  F_mat );
        
       #ifdef VERBOSE  
           std::cout << "GetValue ( " ;
           std::cout << position ;
           std::cout << ") = ";
           std::cout << val.size() << "\n";
           std::cout << val << "\n\n"; 
       #endif    
        
       return val.get(0,0);
    }
    
    /*
     *  Derivative of the objective function is
     *   
     *   \nabla f = A * position - b
     *
     *
     */
    
    void GetDerivative( const ParametersType & position, DerivativeType & derivative ) const
    {
        
      #ifdef VERBOSE
         std::cout << "GetDerivative ( " ;
         for(int j=0; j < m_Dimension-1; ++j) std::cout << position[j] <<" , ";
         std::cout << position[m_Dimension -1] << ") "<< "\n\n";
      #endif    
     
       derivative = DerivativeType(m_Dimension);
       
       //position column vector 
       MatrixType vec(m_Dimension,1);
       for( int i = 0; i < m_Dimension; ++i) vec.put(i,1, position[i]);
       
       MatrixType temp = A * vec -  F_mat;
       for( int i = 0; i < m_Dimension; ++i) derivative[i] = temp.get(i,1);
              
      #ifdef VERBOSE
        std::cout << "(" ; 
        for(int j=0; j < m_Dimension-1; ++j) std::cout << derivative[j] <<" , ";
        std::cout << derivative[m_Dimension -1] << ")" << std::endl;
      #endif     
    }
    
    unsigned int GetNumberOfParameters(void) const
    {
        return m_Dimension;
    }
    
private:
    
    MatrixType A;   //matrix representing the sampled basis
    VectorType F;   // function values
    MatrixType F_mat; 
    
    int m_Dimension;   //number of parameters
    
};//end conjugateCostFunction class



class CommandIterationUpdateConjugateGradient : public itk::Command 
{
public:
  typedef  CommandIterationUpdateConjugateGradient   Self;
  typedef  itk::Command                              Superclass;
  typedef  itk::SmartPointer<Self>                   Pointer;
  itkNewMacro( Self );
  
protected:
  CommandIterationUpdateConjugateGradient() 
  {
    m_IterationNumber=0;
  }
  
public:
  typedef itk::ConjugateGradientOptimizer            OptimizerType;
  typedef const OptimizerType                       *OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( object );
      if( m_FunctionEvent.CheckEvent( &event ) )
        {
        std::cout << m_IterationNumber++ << "   ";
        std::cout << optimizer->GetCachedValue() << "   ";
        std::cout << optimizer->GetCachedCurrentPosition() << std::endl;
        }
      else if( m_GradientEvent.CheckEvent( &event ) )
        {
        std::cout << "Gradient " << optimizer->GetCachedDerivative() << "   ";
        }
    }
    
private:
  unsigned long m_IterationNumber;

  itk::FunctionEvaluationIterationEvent m_FunctionEvent;
  itk::GradientEvaluationIterationEvent m_GradientEvent;
}; //end CommandIterationUpdateConjugateGradient 

namespace itk
{

namespace Functor
{

template< class ImageType, typename TComponentType >
class RshInterpolateFunction
{

public: 

        itkStaticConstMacro(ImageDimension, unsigned int, ImageType::Dimension );
        
        typedef ImageType::PixelType                                                          RshPixelType;        
                                                                                                                                                                                     
        typedef itk::Vector<TComponentType, ImageDimension>                                   PointType;     
        typedef ImageType::IndexType                                                          IndexType;
        
        typedef itk::NeighborhoodIterator<ImageType>                                          NeighborhoodIterator;
        typedef itk::ConjugateGradientOptimizer                                               OptimizerType;
       
        typedef conjugateCostFunction::VectorType                                             VectorType;
        typedef conjugateCostFunction::MatrixType                                             MatrixType;
       
        // Declaration of the CostFunction adaptor
         
	 
/***Constructor****/ 
RshInterpolateFunction(): F_Tolerance(1e-3), G_Tolerance(1e-4), X_Tolerance(1e-8), Epsilon_Function(1e-10), Max_Iterations(100)   
{
  
    m_CostFunction = conjugateCostFunction::New();
    m_Optimizer = OptimizerType::New();
    
    m_Optimizer->SetCostFunction( costFunction.GetPointer() );
    
    //A radius of 1 in all axial directions gives a 3x3x3x3x... neighborhood.
    for (unsigned int i = 0; i < ImageDimension; ++i) radius[i] = 3;  
    m_count = 0; 	    
}     

/***Set the InputImage and the Point, in the image space for interpolation***/
void SetInputImage(const ImageType &in)
{
     m_InputImage = in;  
}	
    
void ComputeCoefficients( typename OptimizerType::Pointer opt) 
{
  typedef  OptimizerType::InternalOptimizerType  vnlOptimizerType;

  // Declaration of the CostFunction adaptor
  conjugateCostFunction::Pointer costFunction = conjugateCostFunction::New();

  opt->SetCostFunction( costFunction.GetPointer() );
  vnlOptimizerType * vnlOp         vnlOptimizer = opt->GetOptimizer();

  vnlOptimizer->set_f_tolerance( F_Tolerance );
  vnlOptimizer->set_g_tolerance( G_Tolerance );
  vnlOptimizer->set_x_tolerance( X_Tolerance ); 
  vnlOptimizer->set_epsilon_function( Epsilon_Function );
  vnlOptimizer->set_max_function_evals( Max_Iterations );

  vnlOptimizer->set_check_derivatives( opt->GetCostFunction()->GetNumberOfParameters() );

  // We start the CG iteration from [100,100,..,100] 
  OptimizerType::ParametersType initialValue (m_Dimension );       // constructor requires vector size
  for(int i=0; i<m_Dimension; ++i) initialValue[i] = 100;

  OptimizerType::ParametersType currentValue(m_Dimension);
  currentValue = initialValue;
  opt->SetInitialPosition( currentValue );

  CommandIterationUpdateConjugateGradient::Pointer observer = CommandIterationUpdateConjugateGradient::New();
  opt->AddObserver( itk::IterationEvent(), observer );
  opt->AddObserver( itk::FunctionEvaluationIterationEvent(), observer );

  try 
    {
    	    opt->StartOptimization();
    }
  catch( itk::ExceptionObject & e )
    {
	    std::cout << "Exception thrown ! " << std::endl;
	    std::cout << "An error ocurred during Optimization" << std::endl;
	    std::cout << "Location    = " << e.GetLocation()    << std::endl;
	    std::cout << "Description = " << e.GetDescription() << std::endl;
	    return EXIT_FAILURE;
    }

#ifdef VERBOSE
  std::cout << "Number of iters = " << opt->GetCurrentIteration()  << std::endl;
  std::cout << "Number of evals = " << vnlOptimizer->get_num_evaluations() << std::endl;    
  std::cout << "Report from vnl optimizer: " << std::endl;
  std::cout << "Stop description   = " << opt->GetStopConditionDescription() << std::endl;
#endif
}//end computeCoefficients

        
/** Do the BSpline Interpolation here
 *  
 *  
 *  points are in image space, not anatomical space
 *
 * 
 */
inline TComponentType operator() ( const PointType & point ) 
{

    RshPixelType tempPixel; 
    bool isInBounds;

    //Find closest discrete point, in the image space
    IndexType discrete_point; 
    for(unsigned int i = 0; i < ImageDimension; ++i) discrete_point[i] = ceil( point[i] );     
    
    NeighborhoodIterator a_it(radius, m_InputImage, m_InputImage->GetRequestedRegion() );
    NeighborhoodIterator b_it(radius, m_InputImage, m_InputImage->GetRequestedRegion() );
    
    a_it.SetPixelPointer( discrete_point );
    b_it.SetPixelPointer( discrete_point );
    
    //fill vnl_vector representing function samples
    std::vector<VectorType::element_type>       vec;
    m_count = 0;
    for (int i = 0; i < a_it.Size(); ++i)
    {
        tempPixel = a_it.GetPixel(i, isInBounds); 
        vec.push_back( tempPixel[3]);
        ++m_count;
    }
   
    VectorType F( vec.size() ); 
    for(int j = 0; j < vec.size(); ++j) F[j] = vec[j];  
     
    //compute BSpline matrix.
    #ifdef VERBOSE
      std::cout<< "vector size: "<< m_count << "\n";
    #endif  
    MatrixType BSplineMat(m_count, m_count);
    ComputeBSplineMat(mat, a_it, b_it);  
        
    //Run CG method on input values
    conjugateCostFunction::Pointer costFunction( BSplineMat, F) = conjugateCostFunction::New();
    m_Optimizer->SetCostFunction( costFunction.GetPointer() );
    ComputeCoefficients( m_Optimizer );
    OptimizerType::ParametersType finalPosition = m_Optimizer->GetCurrentPosition();  //get coefficients
    
    //Value of interpolation evaluated at 'point'.
    TComponentType interpolValue  =  BSplineFunctionInterpolation( point, finalPosition, a_it );
    
    // Get the final value of the optimizer
    #ifdef VERBOSE  
        std::cout << "Testing GetValue() : ";
        OptimizerType::MeasureType finalValue = m_Optimizer->GetValue();
    #endif       
    
    return static_cast<TComponentType>( interpolValue );
}
        
private:

	const double F_Tolerance;  // Function value tolerance
	const double G_Tolerance;  // Gradient magnitude tolerance 
	const double X_Tolerance;  // Search space tolerance
	const double Epsilon_Function; // Step
	const int    Max_Iterations; // Maximum number of iterations
//	const int     SplineOrder;   // Order of Spline function
	
	typename typedef conjugateCostFunction::Pointer                                        m_CostFunction;
	typename typedef OptimizerType::Pointer                                                m_Optimizer;
	
	
	ImageType m_InputImage;
	PointType m_Point; 
	InputImage::IndexType index; 
        NeighborhoodIterator::RadiusType  radius;
        
        int m_count; 
       
    //BSpline function.
    TComponentType BSplineBasisFunction(const TComponentType x, const int k, const TComponentType begin, const TComponentType end);
    
    //Computes the BSpline Matrix 
    void ComputeBSplineMat(MatrixType &mat, const NeighborhoodIterator a, const NeighborhoodIterator b );
   
    //Computes the BSpline interpolation given  
    TComponentType BSplineFunctionInterpolation(PointType &point, OptimizerType::ParametersType finalPosition, const NeighborhoodIterator a);
    
};//end class

} //end namespace Functor

}//end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRshInterpolateImageFunction.txx"
#endif

#endif

