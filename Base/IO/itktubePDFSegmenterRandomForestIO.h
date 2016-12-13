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

#ifndef __itktubePDFSegmenterRandomForestIO_h
#define __itktubePDFSegmenterRandomForestIO_h

#include "itktubeMetaClassPDF.h"

#include "itktubePDFSegmenterRandomForest.h"

namespace itk
{

namespace tube
{

/**
*
* Reads and Writes PDFSegmenterRandomForestIO Files, typically designated .mnda
* files
*
* \author Stephen R. Aylward
*
* \date August 29, 2013
*
*/
template< class TImage, class TLabelMap >
class PDFSegmenterRandomForestIO
{
public:

  typedef PDFSegmenterRandomForestIO<TImage, TLabelMap>  PDFSegmenterIOType;

  typedef PDFSegmenterRandomForest<TImage, TLabelMap>    PDFSegmenterType;

  PDFSegmenterRandomForestIO( void );

  PDFSegmenterRandomForestIO( const char * _headerName );

  PDFSegmenterRandomForestIO( const typename
    PDFSegmenterType::Pointer & _filter );

  ~PDFSegmenterRandomForestIO( void );

  virtual void PrintInfo( void ) const;

  virtual void CopyInfo( const PDFSegmenterIOType & _filterIO );

  virtual void Clear( void );

  virtual bool InitializeEssential( const typename
    PDFSegmenterType::Pointer & _filter );

  void SetPDFSegmenter( const typename
    PDFSegmenterType::Pointer & _filter );

  const typename PDFSegmenterType::Pointer GetPDFSegmenter( void ) const;

  virtual bool CanRead( const char * _headerName = NULL ) const;

  virtual bool Read( const char * _headerName = NULL );

  virtual bool Write( const char * _headerName = NULL );

protected:

  typename PDFSegmenterType::Pointer  m_PDFSegmenter;

}; // End class PDFSegmenterRandomForestIO

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubePDFSegmenterRandomForestIO.hxx"
#endif

#endif // End !defined( __itktubePDFSegmenterRandomForestIO_h )
