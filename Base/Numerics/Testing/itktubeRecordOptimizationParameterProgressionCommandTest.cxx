/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itktubeRecordOptimizationParameterProgressionCommand.h"

#include <itkRegularStepGradientDescentOptimizer.h>

// Taken from itkRegularStepGradientDescentOptimizerTest.cxx
/**
 *  The objectif function is the quadratic form:
 *
 *  1/2 x^T A x - b^T x
 *
 *  Where A is a matrix and b is a vector
 *  The system in this example is:
 *
 *     | 3  2 ||x|   | 2|   |0|
 *     | 2  6 ||y| + |-8| = |0|
 *
 *
 *   the solution is the vector | 2 -2 |
 *
 * \class RSGCostFunction
 */
class RSGCostFunction: public itk::SingleValuedCostFunction
{
public:
  typedef RSGCostFunction               Self;
  typedef itk::SingleValuedCostFunction Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  itkNewMacro( Self );

  enum { SpaceDimension = 2 };

  typedef Superclass::ParametersType      ParametersType;
  typedef Superclass::DerivativeType      DerivativeType;
  typedef Superclass::MeasureType         MeasureType;

  RSGCostFunction()
  {
  }

  MeasureType  GetValue( const ParametersType & parameters ) const
  {
    const double x = parameters[0];
    const double y = parameters[1];

    std::cout << "GetValue( ";
    std::cout << x << " ";
    std::cout << y << " ) = ";

    const MeasureType measure = 0.5*( 3*x*x+4*x*y+6*y*y ) - 2*x + 8*y;

    std::cout << measure << std::endl;
    return measure;
  }

  void GetDerivative( const ParametersType & parameters,
                            DerivativeType  & derivative ) const
  {
    const double x = parameters[0];
    const double y = parameters[1];

    std::cout << "GetDerivative( ";
    std::cout << x << " ";
    std::cout << y << " ) = ";

    derivative = DerivativeType( SpaceDimension );
    derivative[0] = 3 * x + 2 * y -2;
    derivative[1] = 2 * x + 6 * y +8;
    std::cout << derivative[0] << " " << derivative[1] << std::endl;
  }

  unsigned int GetNumberOfParameters( void ) const
    {
    return SpaceDimension;
    }
};

int itktubeRecordOptimizationParameterProgressionCommandTest( int argc, char * argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0]
              << " optimizationParameterProgressionFile" << std::endl;
    return EXIT_FAILURE;
    }
  const char * optimizationParameterProgressionFile = argv[1];

  typedef RSGCostFunction CostFunctionType;
  CostFunctionType::Pointer costFunction = CostFunctionType::New();

  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetCostFunction( costFunction );

  const unsigned int spaceDimension =
                      costFunction->GetNumberOfParameters();

  // We start not so far from  | 2 -2 |
  typedef CostFunctionType::ParametersType ParametersType;
  ParametersType  initialPosition( spaceDimension );
  initialPosition[0] =  100;
  initialPosition[1] = -100;

  typedef OptimizerType::ScalesType ScalesType;
  ScalesType    parametersScale( spaceDimension );
  parametersScale[0] = 1.0;
  parametersScale[1] = 1.0;

  optimizer->MinimizeOn();
  optimizer->SetScales( parametersScale );
  optimizer->SetGradientMagnitudeTolerance( 1e-6 );
  optimizer->SetMaximumStepLength( 30.0 );
  optimizer->SetMinimumStepLength( 1e-6 );
  optimizer->SetNumberOfIterations( 900 );

  optimizer->SetInitialPosition( initialPosition );

  // Record the optimization parameter progression and write to a file.
  const unsigned int NumberOfParameters = CostFunctionType::SpaceDimension;
  typedef itk::tube::RecordOptimizationParameterProgressionCommand< NumberOfParameters >
    RecordOptimizationParameterProgressionCommandType;
  RecordOptimizationParameterProgressionCommandType::Pointer
    recordOptimizationParameterProgressionCommand =
    RecordOptimizationParameterProgressionCommandType::New();
  recordOptimizationParameterProgressionCommand->SetFileName( 
    optimizationParameterProgressionFile );
  recordOptimizationParameterProgressionCommand->Observe( optimizer );
  // some fake fixed parameters
  RecordOptimizationParameterProgressionCommandType::FixedParametersType fixedParameters( 3 );
  fixedParameters[0] = 33.45;
  fixedParameters[1] = 123.4;
  fixedParameters[2] = 746.2;
  recordOptimizationParameterProgressionCommand->SetFixedParameters( fixedParameters );

  try
    {
    optimizer->StartOptimization();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Exception thrown ! " << std::endl;
    std::cerr << "An error occurred during Optimization" << std::endl;
    std::cerr << "Location    = " << e.GetLocation()    << std::endl;
    std::cerr << "Description = " << e.GetDescription() << std::endl;
    return EXIT_FAILURE;
    }

  ParametersType finalPosition = optimizer->GetCurrentPosition();
  std::cout << "Solution        = ( ";
  std::cout << finalPosition[0] << ",";
  std::cout << finalPosition[1] << " )" << std::endl;

  // Make the file was written correctly and has content.
  H5::H5File h5file( optimizationParameterProgressionFile, H5F_ACC_RDONLY );
  H5::DataSet dataset = h5file.openDataSet( "OptimizationParameterProgression" );
  H5::DataSpace dataspace = dataset.getSpace();
  hsize_t numberOfIterations[1];
  dataspace.getSimpleExtentDims( numberOfIterations, NULL );
  std::cout << "Iterations saved: " << numberOfIterations[0] << std::endl;
  const unsigned int expectedNumberOfIterations = 36;
  if( numberOfIterations[0] != expectedNumberOfIterations )
    {
    std::cerr << "The saved iterations, " << numberOfIterations[0]
              << ", does not equal the expected value: "
              << expectedNumberOfIterations << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
