/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include "itktubeSingleValuedCostFunctionImageSource.h"

#include <itkSingleValuedCostFunction.h>
#include <itkCastImageFilter.h>
#include <itkImageFileWriter.h>

namespace itk
{

/** 3D cost function that fades as it travels away from its origin. */
class ArthurDentCostFunction:
  public SingleValuedCostFunction
{
public:
  /** Standard class typedefs. */
  typedef ArthurDentCostFunction        Self;
  typedef SingleValuedCostFunction      Superclass;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro( ArthurDentCostFunction);

  itkNewMacro( Self );

  typedef Superclass::MeasureType         MeasureType;
  typedef Superclass::ParametersType      ParametersType;
  typedef Superclass::ParametersValueType ParametersValueType;

  virtual MeasureType GetValue( const ParametersType & parameters ) const
    override
    {
    return 2 * parameters[0] * parameters[0] + \
           3 * parameters[1] * parameters[1] + \
           4 * parameters[2] * parameters[2];
    }

  virtual void GetDerivative( const ParametersType &,
                             DerivativeType & ) const override
    {
    // irrevelant for this test
    }

  virtual unsigned int GetNumberOfParameters( void ) const override
    {
    return 3;
    }

protected:
  ArthurDentCostFunction() {}
  virtual ~ArthurDentCostFunction() {}

private:
  ArthurDentCostFunction( const Self & ); //purposely not implemented
  void operator=( const Self & );         //purposely not implemented
};

} // end namespace itk

int itktubeSingleValuedCostFunctionImageSourceTest( int argc, char * argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Missing arguments: "
              << argv[0]
              << " outputCostFunctionImage"
              << std::endl;
    return EXIT_FAILURE;
    }
  const char * outputCostFunctionImage = argv[1];

  typedef float      FloatType;
  const unsigned int Dimension = 3;

  typedef itk::ArthurDentCostFunction CostFunctionType;
  CostFunctionType::Pointer costFunction = CostFunctionType::New();

  typedef itk::tube::SingleValuedCostFunctionImageSource< CostFunctionType, Dimension >
    CostFunctionImageSourceType;
  CostFunctionImageSourceType::Pointer costFunctionImageSource =
    CostFunctionImageSourceType::New();
  costFunctionImageSource->SetCostFunction( costFunction );

  typedef CostFunctionImageSourceType::ParametersType ParametersType;

  ParametersType parametersLowerBound( Dimension );
  parametersLowerBound[0] = -5.0;
  parametersLowerBound[1] = -5.0;
  parametersLowerBound[2] = -5.0;
  costFunctionImageSource->SetParametersLowerBound( parametersLowerBound );

  ParametersType parametersUpperBound( Dimension );
  parametersUpperBound[0] = 5.0;
  parametersUpperBound[1] = 5.0;
  parametersUpperBound[2] = 5.0;
  costFunctionImageSource->SetParametersUpperBound( parametersUpperBound );

  ParametersType parametersStep( Dimension );
  parametersStep[0] = 0.5;
  parametersStep[1] = 0.5;
  parametersStep[2] = 0.5;
  costFunctionImageSource->SetParametersStep( parametersStep );

  typedef itk::Image< FloatType, Dimension > OutputImageType;

  typedef itk::CastImageFilter< CostFunctionImageSourceType::OutputImageType,
    OutputImageType > CastImageFilterType;
  CastImageFilterType::Pointer castImageFilter = CastImageFilterType::New();
  castImageFilter->SetInput( costFunctionImageSource->GetOutput() );

  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputCostFunctionImage );
  writer->SetInput( castImageFilter->GetOutput() );
  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught while initializing metric." << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
