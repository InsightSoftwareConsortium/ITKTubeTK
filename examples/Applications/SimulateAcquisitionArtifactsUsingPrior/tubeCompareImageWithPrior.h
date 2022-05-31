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

#ifndef __tubeCompareImageWithPrior_h
#define __tubeCompareImageWithPrior_h

#include "../CLI/tubeCLIProgressReporter.h"

#include <itkImage.h>
#include <itkRigidImageToImageRegistrationMethod.h>
#include <itkTimeProbesCollectorBase.h>

namespace tube
{

template< class TPixel, unsigned int VDimension >
class CompareImageWithPrior
{
public:

  typedef float                                 PixelType;
  typedef itk::Image< PixelType, VDimension >   ImageType;

  typedef itk::RigidImageToImageRegistrationMethod< ImageType >
                                                    RegistrationMethodType;

  CompareImageWithPrior( void );
  ~CompareImageWithPrior( void );

  void SetInput( typename ImageType::Pointer volImage );
  typename ImageType::Pointer GetInput( void );

  void SetMaskImage( typename ImageType::Pointer maskImage );
  typename ImageType::Pointer GetMaskImage( void );

  typename ImageType::Pointer GetOutput( void );
  typename ImageType::Pointer GetOutputMaskImage( void );

  void SetMetricMask( typename ImageType::Pointer metricMask );
  typename ImageType::Pointer GetMetricMask( void );

  void SetForeground( float foreground );

  void SetBackground( float background );

  void SetErode( int erode );

  void SetDilate( int dilate );

  void SetGaussianBlur( float gaussianBlur );

  void SetNormalize( bool normalize );

  void SetSamplingRate( float samplingRate );

  void SetSeed( unsigned int seed );

  void SetUseRegistration( bool reg );
  bool GetUseRegistration( void );
  void SetUseRegistrationTransform( bool reg );
  bool GetUseRegistrationTransform( void );
  void SetUseRegistrationOptimization( bool reg );
  bool GetUseRegistrationOptimization( void );
  void SetRegistrationTransform(
    typename RegistrationMethodType::TransformType::Pointer tfm );
  typename RegistrationMethodType::TransformType::Pointer
    GetRegistrationTransform( void );

  void SetBoundarySize( std::vector< int > & boundarySize );

  void SetTimeCollector( itk::TimeProbesCollectorBase * timeCollector );

  void SetProgressReporter( CLIProgressReporter * progessReporter,
                            float progressStart,
                            float progressRange );

  void Update( void );

  float GetGoodnessOfFit( void );

private:

  typename ImageType::Pointer m_VolImage;
  typename ImageType::Pointer m_MaskImage;
  typename ImageType::Pointer m_OutputVolImage;
  typename ImageType::Pointer m_OutputMaskImage;

  typename ImageType::Pointer m_MetricMask;
  float                       m_Foreground;
  float                       m_Background;
  int                         m_Erode;
  int                         m_Dilate;
  float                       m_GaussianBlur;
  bool                        m_UseRegistration;
  bool                        m_UseRegistrationOptimization;
  bool                        m_UseRegistrationTransform;
  typename RegistrationMethodType::TransformType::Pointer
                              m_RegistrationTransform;
  bool                        m_Normalize;
  std::vector< int >          m_BoundarySize;
  float                       m_SamplingRate;
  unsigned int                m_Seed;

  itk::TimeProbesCollectorBase * m_TimeCollector;
  CLIProgressReporter          * m_ProgressReporter;
  float                          m_ProgressStart;
  float                          m_ProgressRange;

  float                          m_GoF;

}; // End class CompareImageWithPrior

} // End namespace tube

#include "tubeCompareImageWithPrior.hxx"

#endif // End !defined( __tubeCompareImageWithPrior_h )
