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

/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itktubeSpatialObjectToSpatialObjectFilter.h,v $
  Language:  C++
  Date:      $Date: 2004/10/13 15:45:11 $
  Version:   $Revision: 1.4 $
  Author:    Julien Jomier

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itktubeSpatialObjectToSpatialObjectFilter_h
#define __itktubeSpatialObjectToSpatialObjectFilter_h

#include <itkConceptChecking.h>
#include <itkProcessObject.h>

namespace itk
{

namespace tube
{

/**
 * Filter that takes a spatial object as input and produces a spatial object as
 * output.
 *
 * \ingroup  ObjectDocuments
 */
template< class TInputSpatialObject, class TOutputSpatialObject >
class SpatialObjectToSpatialObjectFilter : public ProcessObject
{
public:

  typedef SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
                                              TOutputSpatialObject >
                                      Self;
  typedef ProcessObject               Superclass;
  typedef SmartPointer< Self >        Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  typedef TInputSpatialObject         InputSpatialObjectType;
  typedef TOutputSpatialObject        OutputSpatialObjectType;
  typedef typename InputSpatialObjectType::Pointer
                                      InputSpatialObjectPointer;
  typedef typename InputSpatialObjectType::ConstPointer
                                      InputSpatialObjectConstPointer;

  itkNewMacro( Self );
  itkTypeMacro( SpatialObjectToSpatialObjectFilter, ProcessObject );

  /** The dimension of the input. */
  itkStaticConstMacro( InputSpatialObjectDimension, unsigned int,
                       InputSpatialObjectType::ObjectDimension );

  /** The dimension of the output. */
  itkStaticConstMacro( OutputSpatialObjectDimension, unsigned int,
                       OutputSpatialObjectType::ObjectDimension );

  /** Return the input. */
  virtual const InputSpatialObjectType * GetInput( void );

  /** Return the specified indexed input. */
  virtual const InputSpatialObjectType * GetInput( unsigned int index );

  /** Set the input. */
  using Superclass::SetInput;
  virtual void SetInput( const InputSpatialObjectType * object );

  /** Set the specified indexed input. */
  virtual void SetInput( unsigned int index,
                         const InputSpatialObjectType * object );

protected:

  /** Constructor. */
  SpatialObjectToSpatialObjectFilter( void );

  /** Destructor. */
  virtual ~SpatialObjectToSpatialObjectFilter( void );

  /** Generate the output data. */
  virtual void GenerateData( void );

  /** Generate the information describing the output data. */
  virtual void GenerateOutputInformation( void );

private:

  // Copy constructor not implemented.
  SpatialObjectToSpatialObjectFilter( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

}; // End class SpatialObjectToSpatialObjectFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSpatialObjectToSpatialObjectFilter.hxx"
#endif

#endif // End !defined(__itktubeSpatialObjectToSpatialObjectFilter_h)
