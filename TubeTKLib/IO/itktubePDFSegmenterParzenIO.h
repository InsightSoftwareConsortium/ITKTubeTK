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

#ifndef __itktubePDFSegmenterParzenIO_h
#define __itktubePDFSegmenterParzenIO_h

#include "itktubeMetaClassPDF.h"

#include "itktubePDFSegmenterParzen.h"

namespace itk
{

namespace tube
{

/**
*
* Reads and Writes PDFSegmenterParzenIO Files, typically designated .mnda
* files
*
* \author Stephen R. Aylward
*
* \date August 29, 2013
*
*/
template< class TImage, class TLabelMap >
class PDFSegmenterParzenIO
{
public:

  typedef PDFSegmenterParzenIO< TImage, TLabelMap >  PDFSegmenterIOType;

  typedef PDFSegmenterParzen< TImage, TLabelMap >    PDFSegmenterType;

  PDFSegmenterParzenIO( void );

  PDFSegmenterParzenIO( const char * _headerName );

  PDFSegmenterParzenIO( const typename
    PDFSegmenterType::Pointer & _filter );

  ~PDFSegmenterParzenIO( void );

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

}; // End class PDFSegmenterParzenIO

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubePDFSegmenterParzenIO.hxx"
#endif

#endif // End !defined( __itktubePDFSegmenterParzenIO_h )
