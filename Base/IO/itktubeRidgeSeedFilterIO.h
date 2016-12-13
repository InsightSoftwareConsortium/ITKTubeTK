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

#ifndef __itktubeRidgeSeedFilterIO_h
#define __itktubeRidgeSeedFilterIO_h

#include "itktubeMetaRidgeSeed.h"

#include "itktubeRidgeSeedFilter.h"

namespace itk
{

namespace tube
{

/**
*
* Reads and Writes RidgeSeedFilterIO Files, typically designated .mnda files
*
* \author Stephen R. Aylward
*
* \date August 29, 2013
*
*/
template< class TImage, class TLabelMap >
class RidgeSeedFilterIO
{
public:

  typedef RidgeSeedFilterIO< TImage, TLabelMap >  RidgeSeedFilterIOType;

  typedef RidgeSeedFilter< TImage, TLabelMap >    RidgeSeedFilterType;

  RidgeSeedFilterIO( void );

  RidgeSeedFilterIO( const char * _headerName );

  RidgeSeedFilterIO( const typename
    RidgeSeedFilterType::Pointer & _filter );

  ~RidgeSeedFilterIO( void );

  virtual void PrintInfo( void ) const;

  virtual void CopyInfo( const RidgeSeedFilterIOType & _filterIO );

  virtual void Clear( void );

  virtual bool InitializeEssential( const typename
    RidgeSeedFilterType::Pointer & _filter );

  void SetRidgeSeedFilter( const typename
    RidgeSeedFilterType::Pointer & _filter );

  const typename RidgeSeedFilterType::Pointer GetRidgeSeedFilter( void ) const;

  virtual bool CanRead( const char * _headerName = NULL ) const;

  virtual bool Read( const char * _headerName = NULL );

  virtual bool Write( const char * _headerName = NULL );

protected:

  typename RidgeSeedFilterType::Pointer  m_RidgeSeedFilter;

}; // End class RidgeSeedFilterIO

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeRidgeSeedFilterIO.hxx"
#endif

#endif // End !defined( __itktubeRidgeSeedFilterIO_h )
