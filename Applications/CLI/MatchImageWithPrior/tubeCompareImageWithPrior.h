/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved. 

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __tubeCompareImageWithPrior_h
#define __tubeCompareImageWithPrior_h

#include "itkOrientedImage.h"
#include "itkRigidImageToImageRegistrationMethod.h"
#include "itkTimeProbesCollectorBase.h"

#include "tubeCLIProgressReporter.h"

namespace tube
{

template< class pixelT, unsigned int dimensionT >
class CompareImageWithPrior
{
public:

  typedef float                                         PixelType;
  typedef itk::OrientedImage< PixelType,  dimensionT >  ImageType;

  typedef itk::RigidImageToImageRegistrationMethod< ImageType >
                                                    RegistrationMethodType;

  CompareImageWithPrior( void );
  ~CompareImageWithPrior( void );

  void SetVolumeImage( typename ImageType::Pointer volImage );
  typename ImageType::Pointer GetVolumeImage( void );

  void SetMaskImage( typename ImageType::Pointer maskImage );
  typename ImageType::Pointer GetMaskImage( void );

  void SetOriginalMaskImage( typename ImageType::Pointer orgMaskImage );
  typename ImageType::Pointer GetOriginalMaskImage( void );

  void SetForeground( float foreground );

  void SetErode( int erode );

  void SetDilate( int dilate );

  void SetGaussianBlur( float gaussianBlur );

  void SetUseRegistration( bool reg );
  void SetUseRegistrationTransform( bool reg );
  void SetRegistrationTransform( 
    typename RegistrationMethodType::TransformType::Pointer tfm );
  typename RegistrationMethodType::TransformType::Pointer 
    GetRegistrationTransform( void );

  void SetBoundarySize( std::vector< int > & boundarySize );

  void SetUseMeanSquaresMetric( bool useMeanSquaresMetric );

  void SetTimeCollector( itk::TimeProbesCollectorBase * timeCollector );

  void SetProgressReporter( CLIProgressReporter * progessReporter,
                            float progressStart,
                            float progressRange );

  void Update( void );

  float GetGoodnessOfFit( void );

private:

  typename ImageType::Pointer m_VolImage;
  typename ImageType::Pointer m_MaskImage;
  typename ImageType::Pointer m_OrgMaskImage;
  float                       m_Foreground;
  int                         m_Erode;
  int                         m_Dilate;
  float                       m_GaussianBlur;
  bool                        m_UseRegistration;
  bool                        m_UseRegistrationTransform;
  typename RegistrationMethodType::TransformType::Pointer
                              m_RegistrationTransform;
  bool                        m_Normalize;
  std::vector< int >          m_BoundarySize;
  bool                        m_UseMeanSquaresMetric;
  bool                        m_UseCorrelationMetric;
  float                       m_SamplingRate;

  itk::TimeProbesCollectorBase * m_TimeCollector;
  CLIProgressReporter          * m_ProgressReporter;
  float                          m_ProgressStart;
  float                          m_ProgressRange;

  float                          m_GoF;

};

}

#include "tubeCompareImageWithPrior.txx"

#endif
