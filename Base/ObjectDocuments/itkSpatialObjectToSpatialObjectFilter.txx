/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkSpatialObjectToSpatialObjectFilter.txx,v $
  Language:  C++
  Date:      $Date: 2004/08/11 20:59:14 $
  Version:   $Revision: 1.2 $
  Author:    Julien Jomier

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkSpatialObjectToSpatialObjectFilter_txx
#define __itkSpatialObjectToSpatialObjectFilter_txx

#include "itkSpatialObjectToSpatialObjectFilter.h"


namespace itk
{

namespace tube
{


template <class TInputSpatialObject, class TOutputSpatialObject>
SpatialObjectToSpatialObjectFilter<TInputSpatialObject,TOutputSpatialObject>
::SpatialObjectToSpatialObjectFilter( void )
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(1);
}


template <class TInputSpatialObject, class TOutputSpatialObject>
SpatialObjectToSpatialObjectFilter<TInputSpatialObject,TOutputSpatialObject>
::~SpatialObjectToSpatialObjectFilter( void )
{
}


template <class TInputSpatialObject, class TOutputSpatialObject>
void
SpatialObjectToSpatialObjectFilter<TInputSpatialObject,TOutputSpatialObject>
::SetInput(const InputSpatialObjectType *input)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(0,
    const_cast< InputSpatialObjectType * >( input ) );
}


/**
 * Connect one of the operands for pixel-wise addition
 */
template <class TInputSpatialObject, class TOutputSpatialObject>
void
SpatialObjectToSpatialObjectFilter<TInputSpatialObject,TOutputSpatialObject>
::SetInput( unsigned int index, const TInputSpatialObject * object )
{
  if( index+1 > this->GetNumberOfInputs() )
    {
    this->SetNumberOfRequiredInputs( index + 1 );
    }

  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(index,
     const_cast< TInputSpatialObject *>( object ) );
}


template <class TInputSpatialObject, class TOutputSpatialObject>
const typename SpatialObjectToSpatialObjectFilter<TInputSpatialObject,TOutputSpatialObject>::InputSpatialObjectType *
SpatialObjectToSpatialObjectFilter<TInputSpatialObject,TOutputSpatialObject>
::GetInput( void )
{
  if(this->GetNumberOfInputs() < 1)
    {
    return 0;
    }

  return static_cast<const TInputSpatialObject * >
    (this->ProcessObject::GetInput(0) );
}


template <class TInputSpatialObject, class TOutputSpatialObject>
const typename SpatialObjectToSpatialObjectFilter<TInputSpatialObject,TOutputSpatialObject>::InputSpatialObjectType *
SpatialObjectToSpatialObjectFilter<TInputSpatialObject,TOutputSpatialObject>
::GetInput(unsigned int idx)
{
  return static_cast< const TInputSpatialObject * >
    (this->ProcessObject::GetInput(idx));
}


template<class TInputSpatialObject, class TOutputSpatialObject>
void
SpatialObjectToSpatialObjectFilter<TInputSpatialObject,TOutputSpatialObject>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itkSpatialObjectToSpatialObjectFilter_txx)
