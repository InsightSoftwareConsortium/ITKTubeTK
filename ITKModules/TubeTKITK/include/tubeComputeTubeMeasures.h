/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __tubeComputeTubeMeasures_h
#define __tubeComputeTubeMeasures_h

// ITK includes
#include <itkObject.h>

// TubeTK includes
#include <itktubeComputeTubeMeasuresFilter.h>
#include "tubeWrappingMacros.h"

namespace tube
{
/** \class ComputeTubeMeasures
 *
 *  \ingroup TubeTKITK
 */

template< class TPixel, unsigned int Dimension >
class ComputeTubeMeasures:
  public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef ComputeTubeMeasures                        Self;
  typedef itk::Object                                SuperClass;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  typedef itk::tube::ComputeTubeMeasuresFilter
  < TPixel, Dimension > FilterType;

  typedef typename FilterType::InputImageType       InputImageType;
  typedef typename FilterType::OutputImageType      OutputImageType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods ). */
  itkTypeMacro( ComputeTubeMeasures, Object );

  /** Set/Get scale */
  tubeWrapSetMacro( Scale, int, Filter );
  tubeWrapGetMacro( Scale, int, Filter );

  /* Set/Get input image */
  tubeWrapSetConstObjectMacro( InputImage, InputImageType, Filter );
  tubeWrapGetConstObjectMacro( InputImage, InputImageType, Filter );

  /* Runs the application */
  tubeWrapUpdateMacro( Filter );

  /** Get output Ridge Image */
  tubeWrapGetObjectMacro( Ridgeness, OutputImageType, Filter );

  /** Get output Round Image */
  tubeWrapGetObjectMacro( Roundness, OutputImageType, Filter );

  /** Get output Curvature Image */
  tubeWrapGetObjectMacro( Curvature, OutputImageType, Filter );

  /** Get output Levelness Image */
  tubeWrapGetObjectMacro( Levelness, OutputImageType, Filter );

protected:
  ComputeTubeMeasures( void );
  ~ComputeTubeMeasures() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itktubeComputeTubeMeasuresFilter parameters **/
  ComputeTubeMeasures( const Self & );
  void operator=( const Self & );

  typename FilterType::Pointer m_Filter;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeComputeTubeMeasures.hxx"
#endif

#endif // End !defined( __tubeComputeTubeMeasures_h )
