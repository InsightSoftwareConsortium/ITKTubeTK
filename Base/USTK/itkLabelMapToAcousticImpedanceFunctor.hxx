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

#ifndef __itkLabelMapToAcousticImpedanceFunctor_hxx
#define __itkLabelMapToAcousticImpedanceFunctor_hxx

#include "itkLabelMapToAcousticImpedanceFunctor.h"

#include <cstddef>

namespace itk
{
namespace Functor
{

template< class TLabelPixel, class TImpedancePixel, class TLookupTable >
LabelMapToAcousticImpedanceFunctor< TLabelPixel, TImpedancePixel, TLookupTable >
::LabelMapToAcousticImpedanceFunctor( void )
  : m_LookupTable( NULL )
{
}


template< class TLabelPixel, class TImpedancePixel, class TLookupTable >
void
LabelMapToAcousticImpedanceFunctor< TLabelPixel, TImpedancePixel, TLookupTable >
::SetLookupTable( const LookupTableType * lookupTable )
{
  this->m_LookupTable = lookupTable;
}


template< class TLabelPixel, class TImpedancePixel, class TLookupTable >
const typename
LabelMapToAcousticImpedanceFunctor< TLabelPixel, TImpedancePixel, TLookupTable >
::LookupTableType *
LabelMapToAcousticImpedanceFunctor< TLabelPixel, TImpedancePixel, TLookupTable >
::GetLookupTable( void ) const
{
  return this->m_LookupTable;
}


template< class TLabelPixel, class TImpedancePixel, class TLookupTable >
bool
LabelMapToAcousticImpedanceFunctor< TLabelPixel, TImpedancePixel, TLookupTable >
::operator!=( const LabelMapToAcousticImpedanceFunctor & other ) const
{
  return this->m_LookupTable != other.m_LookupTable;
}


template< class TLabelPixel, class TImpedancePixel, class TLookupTable >
bool
LabelMapToAcousticImpedanceFunctor< TLabelPixel, TImpedancePixel, TLookupTable >
::operator==( const LabelMapToAcousticImpedanceFunctor & other ) const
{
    return !( *this != other );
}

} // End namespace Functor

} // End namespace itk

#endif // End !defined( __itkLabelMapToAcousticImpedanceFunctor_hxx )
