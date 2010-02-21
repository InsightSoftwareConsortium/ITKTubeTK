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

#ifndef __itkGeometryObject_txx
#define __itkGeometryObject_txx

#include "itkGeometryObject.h"
#include "itkNumericTraits.h"
#include <algorithm>
#include <string>

namespace itk 
{

namespace Geometry
{

/** Constructor */
template< unsigned int TDimension >
GeometryObject< TDimension >
::GeometryObject( void )
{
}

/** Destructor */
template< unsigned int TDimension >
GeometryObject< TDimension >
::~GeometryObject( void )
{
}


template< unsigned int TDimension >
void
GeometryObject< TDimension >
::CopyInformation(const DataObject * itkNotUsed(data))
{
}

template< unsigned int TDimension >
void
GeometryObject< TDimension >
::UpdateOutputInformation()
{
}

template< unsigned int TDimension >
void
GeometryObject< TDimension >
::Update()
{
}

template< unsigned int TDimension >
void
GeometryObject< TDimension >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // end of namespace geometry

} // end of namespace itk

#endif // __GeometryObject_txx
