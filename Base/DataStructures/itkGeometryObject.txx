/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGeometryObject.txx,v $
  Language:  C++
  Date:      $Date: 2009-01-28 20:10:27 $
  Version:   $Revision: 1.77 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

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
::CopyInformation(const DataObject *data)
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
