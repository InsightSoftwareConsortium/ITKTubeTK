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
#ifndef __tubeSubImageGenerator_h
#define __tubeSubImageGenerator_h

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
class SubImageGenerator
{
public:

  typedef itk::OrientedImage<pixelT,dimensionT> ImageType;

  /// Default Constructor
  SubImageGenerator();

  /// Full populated constructor
  SubImageGenerator( std::vector<int> roiCenter,
                     std::vector<int> roiSize,
                     typename ImageType::Pointer inputVolume,
                     typename ImageType::Pointer inputMask );

  /// Default Destructor
  virtual ~SubImageGenerator();

  /// Update function for doing the processing and producing the output 
  /// volumes.
  void Update();

  void SetRoiCenter( std::vector<int> roiCenter );
  void SetRoiSize( std::vector<int> roiSize );
  void SetInputVolume( typename ImageType::Pointer inputVolume );
  void SetInputMask( typename ImageType::Pointer inputMask );

  typename ImageType::Pointer GetOutputVolume();
  typename ImageType::Pointer GetOutputMask();

protected:
  
  std::vector<int>              m_RoiCenter;
  std::vector<int>              m_RoiSize;
  typename ImageType::Pointer   m_InputVolume;
  typename ImageType::Pointer   m_InputMask;
  typename ImageType::Pointer   m_OutputVolume;
  typename ImageType::Pointer   m_OutputMask;
  
};

}

#include "tubeSubImageGenerator.txx"

#endif
