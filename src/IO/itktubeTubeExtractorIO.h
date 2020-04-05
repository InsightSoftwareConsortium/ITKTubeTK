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

#ifndef __itktubeTubeExtractorIO_h
#define __itktubeTubeExtractorIO_h

#include "itktubeTubeExtractor.h"

#include "itktubeMetaTubeExtractor.h"

namespace itk
{

namespace tube
{

/**
*
* Reads and Writes TubeExtractorIO Files, typically designated .mnda files
*
* \author Stephen R. Aylward
*
* \date August 29, 2013
*
*/
template< class TImage >
class TubeExtractorIO
{
public:

  typedef TubeExtractor< TImage >    TubeExtractorType;

  TubeExtractorIO( void );

  TubeExtractorIO( const char * _headerName );

  TubeExtractorIO( const typename
    TubeExtractorType::Pointer & _filter );

  ~TubeExtractorIO( void );

  virtual void PrintInfo( void ) const;

  virtual void CopyInfo( const TubeExtractorIO< TImage > & _filterIO );

  virtual void Clear( void );

  virtual bool InitializeEssential( const typename
    TubeExtractorType::Pointer & _filter );

  void SetTubeExtractor( TubeExtractorType * _filter );

  const typename TubeExtractorType::Pointer GetTubeExtractor( void ) const;

  virtual bool CanRead( const char * _headerName = NULL ) const;

  virtual bool Read( const char * _headerName = NULL );

  virtual bool Write( const char * _headerName = NULL );

protected:

  typename TubeExtractorType::Pointer  m_TubeExtractor;

}; // End class TubeExtractorIO

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTubeExtractorIO.hxx"
#endif

#endif // End !defined( __itktubeTubeExtractorIO_h )
