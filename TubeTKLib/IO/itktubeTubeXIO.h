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
#ifndef __itktubeTubeXIO_h
#define __itktubeTubeXIO_h

#include "itkVesselTubeSpatialObject.h"
#include "itkGroupSpatialObject.h"

#include <list>

namespace itk
{

namespace tube
{

template< unsigned int TDimension = 3 >
class TubeXIO : public Object
{
public:

  typedef TubeXIO                                 Self;
  typedef Object                                  Superclass;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;

  typedef VesselTubeSpatialObject< TDimension >   TubeType;
  typedef GroupSpatialObject< TDimension >        TubeGroupType;

  typedef Size< TDimension >                      SizeType;

  itkTypeMacro( TubeXIO, Object );

  itkNewMacro( TubeXIO );

  bool  Read( const std::string & _filename );

  bool  Write( const std::string & _filename );

  void  SetTubeGroup( TubeGroupType * _tubes );

  typename TubeGroupType::Pointer & GetTubeGroup( void );

  /** Set the TubeX file dimensions */
  itkSetMacro( Dimensions, SizeType );
  itkGetConstMacro( Dimensions, SizeType );

protected:

  TubeXIO( void );
  virtual ~TubeXIO( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  TubeXIO( const Self& );
  void operator=( const Self& );

  typename TubeGroupType::Pointer  m_TubeGroup;
  SizeType                         m_Dimensions;

}; // TubeXIO

} // namespace tube

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTubeXIO.hxx"
#endif

#endif
