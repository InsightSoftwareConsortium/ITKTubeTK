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
#ifndef __tubeJointHistogramGenerator_h
#define __tubeJointHistogramGenerator_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <vector>

// It is important to use OrientedImages
#include "itkOrientedImage.h"

namespace tube
{

// Function for Getting a good region
template< class pixelT, unsigned int dimensionT >
class JointHistogramGenerator
{
public:

  typedef itk::Image<pixelT,dimensionT> ImageType;

  /// Default Constructor
  JointHistogramGenerator();

  /// Default Destructor
  virtual ~JointHistogramGenerator();

  /// Update function for doing the processing and producing the output 
  /// volumes.
  void Update();

  void SetInputVolume( typename ImageType::Pointer inputVolume );
  void SetInputMask( typename ImageType::Pointer inputMask );

  void SetNumberOfBins( unsigned int numBins );

  typename ImageType::Pointer GetOutputVolume();
  typename ImageType::Pointer GetOutputMask();

protected:
  
  typename ImageType::Pointer   m_InputVolume;
  typename ImageType::Pointer   m_InputMask;
  typename ImageType::Pointer   m_OutputVolume;
  unsigned int                  m_NumberOfBins;
  
};

}

#include "tubeJointHistogramGenerator.txx"

#endif
