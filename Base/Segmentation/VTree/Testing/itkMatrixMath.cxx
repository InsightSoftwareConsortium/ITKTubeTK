#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"

#include "itkMersenneTwisterRandomVariateGenerator.h"

#include "../itkMatrixMath.h"

template< int dimensionT >
int Test( void )
{
  double epsilon = 0.000001;

  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer rndGen
    = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
  rndGen->Initialize( 1 );

  int returnStatus = EXIT_SUCCESS;

  for( unsigned int count=0; count<1000; count++ )
    {
    vnl_vector<float> v1(dimensionT);
    for( unsigned int d=0; d<dimensionT; d++ )
      {
      v1[d] = rndGen->GetVariateWithOpenRange();
      }
    vnl_vector<float> v2(dimensionT);
    if( dimensionT == 3 )
      {
      v2 = itk::GetOrthogonalVector( v1 );
      if( dot_product( v1, v2 ) > epsilon )
        {
        std::cerr << count << " : ";
        std::cerr << "FAILURE: GetOrthogonalVector: DotProduct = "
                  << v1 << " .* " << v2 << " = " 
                  << dot_product( v1, v2 ) << std::endl;
        returnStatus = EXIT_FAILURE;
        }
    
      vnl_vector<float> v3 = itk::GetCrossVector( v1, v2 );
      if( dot_product( v1, v3 ) > epsilon ||
          dot_product( v2, v3 ) > epsilon )
        {
        std::cerr << count << " : ";
        std::cerr << "FAILURE: GetCrossVector: DotProduct = "
          << dot_product( v1, v3 ) << " and " 
          << dot_product( v2, v3 ) << std::endl;
        returnStatus = EXIT_FAILURE;
        }
      }
    else
      {
      for( unsigned int d=0; d<dimensionT; d++ )
        {
        v2[d] = rndGen->GetVariateWithOpenRange();
        }
      }
  
    v2 = v2.normalize();
    vnl_vector<float> v4 = itk::ComputeLineStep( v1, 0.5, v2 );
    if( vnl_math_abs( itk::ComputeEuclideanDistanceVector( v1, v4 ) - 0.5 )
        > epsilon )
      {
      std::cerr << count << " : ";
      std::cerr << "FAILURE: ComputeLineStep = "
        << v1 << " + 0.5 * " << v2 << " != " << v4 << std::endl;
      std::cerr << "FAILURE: ComputeEuclidenDistanceVector = "
        << itk::ComputeEuclideanDistanceVector( v1, v4 ) << std::endl;
      returnStatus = EXIT_FAILURE;
      }

    vnl_matrix<float> m1(dimensionT, dimensionT);
    for( unsigned int r=0; r<dimensionT; r++ )
      {
      for( unsigned int c=r; c<dimensionT; c++ )
        {
        m1(r,c) = rndGen->GetVariateWithOpenRange();
        m1(r,c) = m1(c,r);
        }
      }

    vnl_matrix<float> eVects(dimensionT, dimensionT);
    vnl_vector<float> eVals(dimensionT);
    itk::Eigen( m1, eVects, eVals, true );
    for( unsigned int d=0; d<dimensionT; d++ )
      {
      v1 = m1 * eVects.get_column(d);
      if( vnl_math_abs( v1.magnitude() - eVals[d] ) > epsilon )
        {
        std::cerr << count << " : ";
        std::cerr << "FAILURE: Eigen : "
          << " v1 * M1 = " << v1 
          << " and v1 norm = " << v1.magnitude()
          << " != " << eVals[d] 
          << std::endl;
        returnStatus = EXIT_FAILURE;
        }
      }
    }
  
  return returnStatus;
}

int itkMatrixMathTest(int argc, char * argv[])
{

  if( Test<2>() == EXIT_FAILURE || 
      Test<3>() == EXIT_FAILURE ||
      Test<4>() == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;

}

