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

#ifndef __itkLabelMapToAcousticImpedanceImageFilter_hxx
#define __itkLabelMapToAcousticImpedanceImageFilter_hxx

#include "itkLabelMapToAcousticImpedanceImageFilter.h"

namespace itk
{

template< class TInputImage, class TOutputImage, class TLookupTable >
void
LabelMapToAcousticImpedanceImageFilter< TInputImage, TOutputImage,
  TLookupTable >
::BeforeThreadedGenerateData( void )
{
  Superclass::BeforeThreadedGenerateData();

  const typename FunctorType::LookupTableType * lookupTable =
    this->GetFunctor().GetLookupTable();
  if( lookupTable == NULL )
    {
    itkExceptionMacro( << "Please set the lookup table \
 for the LabelMapToAcousticImpedanceImageFilter functor. " );
    }
}

} // End namespace itk

#endif // End !defined( __itkLabelMapToAcousticImpedanceImageFilter_hxx )
