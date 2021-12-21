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

#ifndef __itktubeEnhanceContrastUsingPriorImageFilter_hxx
#define __itktubeEnhanceContrastUsingPriorImageFilter_hxx


namespace itk
{

namespace tube
{

//----------------------------------------------------------------------------
template< class TPixel, unsigned int VDimension >
EnhanceContrastUsingPriorImageFilter< TPixel, VDimension >
::EnhanceContrastUsingPriorImageFilter( void )
{
  m_InputMaskImage = NULL;
}

//----------------------------------------------------------------------------
template< class TPixel, unsigned int VDimension >
void
EnhanceContrastUsingPriorImageFilter< TPixel, VDimension >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

}

//----------------------------------------------------------------------------
template< class TPixel, unsigned int VDimension >
void
EnhanceContrastUsingPriorImageFilter< TPixel, VDimension >
::GenerateData( void )
{
  const ImageType *inputImage = this->GetInput();
  ImageType *outputImage = this->GetOutput( 0 );

  outputImage->CopyInformation( inputImage );
  outputImage->SetRegions( inputImage->GetLargestPossibleRegion() );
  outputImage->Allocate();

  itk::ImageRegionConstIterator< ImageType > iter( inputImage,
    inputImage->GetLargestPossibleRegion() );
  double inputMin = iter.Get();
  double inputMax = inputMin;
  while( !iter.IsAtEnd() )
    {
    double tf = iter.Get();
    if( tf < inputMin )
      {
      inputMin = tf;
      }
    else if( tf > inputMax )
      {
      inputMax = tf;
      }
    ++iter;
    }

  itk::Array< double > params( 3 );
  params[0] = m_ObjectScale;
  params[1] = m_BackgroundScale;
  params[2] = 20 * ( inputMax - inputMin );

  typename ContrastCostFunctionType::Pointer costFunc =
    ContrastCostFunctionType::New();
  costFunc->SetInputImage( inputImage );
  costFunc->SetInputMask( m_InputMaskImage );
  costFunc->SetOutputImage( outputImage );
  costFunc->SetMaskObjectValue( m_MaskObjectValue );
  costFunc->SetMaskBackgroundValue( m_MaskBackgroundValue );

  InitialOptimizerType::Pointer initOptimizer =
    InitialOptimizerType::New();

  typename itk::Statistics::NormalVariateGenerator::Pointer normGen =
    itk::Statistics::NormalVariateGenerator::New();
  if( m_OptimizationSeed > 0 )
    {
    normGen->Initialize( m_OptimizationSeed );
    }
  initOptimizer->SetNormalVariateGenerator( normGen );
  initOptimizer->Initialize( 1.0 );
  initOptimizer->SetMetricWorstPossibleValue( 101 );
  initOptimizer->SetMaximumIteration( m_OptimizationIterations * 0.5 );
  initOptimizer->SetMaximize( true );

  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetUseUnitLengthGradient( true );
  optimizer->SetMaximumIteration( m_OptimizationIterations * 0.4 );
  optimizer->SetMaximumLineIteration( m_OptimizationIterations * 0.2 );
  optimizer->SetStepLength( 0.1 );
  optimizer->SetStepTolerance( 0.001 );
  optimizer->SetValueTolerance( 0.01 );
  optimizer->SetMaximize( true );

  OptimizerType::ScalesType scales( 3 );
  scales[0] = 1.0 / 0.1;
  scales[1] = 1.0 / 2.0;
  scales[2] = 1.0 / ( params[2] / 10 );

  typename ContrastCostFunctionType::ParametersType costFunctionScales( 3 );
  costFunctionScales[0] = scales[0];
  costFunctionScales[1] = scales[1];
  costFunctionScales[2] = scales[2];

  OptimizerType::ScalesType scales2( 3 );
  scales2[0] = scales[0] * scales[0];
  scales2[1] = scales[1] * scales[1];
  scales2[2] = scales[2] * scales[2];

  // OnePlusOne should be passed squared-scales
  initOptimizer->SetScales( scales2 );
  optimizer->SetScales( scales );
  costFunc->SetScales( costFunctionScales );

  initOptimizer->SetCostFunction( costFunc );
  optimizer->SetCostFunction( costFunc );
  costFunc->SetOutputImage( outputImage );
  costFunc->Initialize();
  initOptimizer->SetInitialPosition( params );

  initOptimizer->StartOptimization();

  params = initOptimizer->GetCurrentPosition();
  double result = costFunc->GetValue( params );
  std::cout << "Intermediate params = " << params
            << " Result = " << result << std::endl;

  optimizer->SetInitialPosition( params );
  optimizer->StartOptimization();

  params = optimizer->GetCurrentPosition();
  result = costFunc->GetValue( params );
  std::cout << "Winning params = " << params
            << " Result = " << result << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeEnhanceContrastUsingPriorImageFilter_hxx )
