/*=========================================================================

Library:   TubeTK

Copyright 2012 Kitware Inc. 28 Corporate Drive,
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
#include <itkLabelMapToAcousticImpedanceImageFilter.h>

#include <map>

int itkLabelMapToAcousticImpedanceImageFilterTest( int argc, char * argv [] )
{
  // Use std::vector here for efficiency if TLabelPixel is always an integer?
  //typedef std::map< TLabelPixel, TImpedancePixel > LookupTableType;
  return EXIT_SUCCESS;
}
