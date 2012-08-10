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
#ifndef __itkTubeMetaLDA_h
#define __itkTubeMetaLDA_h

#include "metaTypes.h"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include "metaUtils.h"
#include "metaForm.h"

/*!    MetaLDA ( .h and .cpp )
*
* Description:
*    Reads and Writes MetaLDA Files, typically designated .mlda files
*
*    REQUIRED: itkTubeLDAGenerator instance
*
* \author Stephen R. Aylward
*
* \date August 29, 1999
*
* Depends on:
*    MetaUtils.h
*    MetaForm.h*/

namespace itk {

namespace tube {

class METAIO_EXPORT MetaLDA
: public MetaForm
{
  /////
  //
  // PUBLIC
  //
  ////
public:

  typedef vnl_vector< double >    LDAValuesType;

  typedef vnl_matrix< double >    LDAMatrixType;

  ////
  //
  // Constructors & Destructor
  //
  ////
  MetaLDA( void );

  MetaLDA( const char *_headerName );

  MetaLDA( const MetaLDA & _metaLDA );

  MetaLDA( const LDAValuesType & _ldaValues,
      const LDAMatrixType & _ldaMatrix );

  ~MetaLDA( void );

  virtual void  PrintInfo( void ) const;

  using MetaForm::CopyInfo;
  virtual void  CopyInfo( const MetaLDA & _lda );

  virtual void  Clear( void );

  bool  InitializeEssential( const LDAValuesType & _ldaValues,
      const LDAMatrixType & _ldaMatrix );

  void  SetLDAValues( const LDAValuesType & _ldaValues );
  const LDAValuesType & GetLDAValues( void ) const;

  void  SetLDAMatrix( const LDAMatrixType & _ldaMatrix );
  const LDAMatrixType & GetLDAMatrix( void ) const;

  virtual bool CanRead( const char *_headerName=NULL ) const;
  virtual bool Read( const char *_headerName=NULL );
  virtual bool CanReadStream( METAIO_STREAM::ifstream * _stream ) const;

  virtual bool ReadStream( METAIO_STREAM::ifstream * _stream );

  virtual bool Write( const char *_headName=NULL );

  virtual bool WriteStream( METAIO_STREAM::ofstream * _stream );

  ////
  //
  // PROTECTED
  //
  ////
protected:

  LDAValuesType   m_LDAValues;

  LDAMatrixType   m_LDAMatrix;

  void  M_Destroy( void );

  void  M_SetupReadFields( void );

  void  M_SetupWriteFields( void );

  bool  M_Read( void );

};

}

}

#endif
