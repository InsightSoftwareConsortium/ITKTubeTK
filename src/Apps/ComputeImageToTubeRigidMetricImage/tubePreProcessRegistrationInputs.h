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

#ifndef __tubePreProcessRegistrationInputs_h
#define __tubePreProcessRegistrationInputs_h

namespace tube {

template<
  unsigned int VDimension,
  typename TFloat,
  typename TTube,
  typename TTubeNet,
  typename TImage,
  typename TRegistrationMethod >
int
PreProcessRegistrationInputs( int argc,
  char * argv[],
  itk::TimeProbesCollectorBase & timeCollector,
  typename TImage::Pointer & currentImage,
  typename TTubeNet::Pointer & tubeNet,
  typename TRegistrationMethod::FeatureWeightsType & pointWeights );

} // namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubePreProcessRegistrationInputs.hxx"
#endif

#endif // End !defined( __tubePreProcessRegistrationInputs_h )
