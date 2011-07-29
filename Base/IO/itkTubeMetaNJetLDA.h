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
#include "metaTypes.h"

#ifndef __tubeMetaNJetLDA_h
#define __tubeMetaNJetLDA_h

#include <vector>

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include "metaUtils.h"
#include "metaForm.h"
#include "itkTubeMetaLDA.h"

/*!    MetaNJetLDA ( .h and .cpp )
*
* Description:
*    Reads and Writes MetaNJetLDA Files, typically designated .mnda files
*
*    REQUIRED: itkTubeNJetLDAGenerator instance
*
* \author Stephen R. Aylward
*
* \date August 29, 1999
*
* Depends on:
*    MetaUtils.h
*    MetaLDA.h
*    MetaForm.h
*/

namespace itk {

namespace tube {

class METAIO_EXPORT MetaNJetLDA
: public MetaLDA
{
/////
//
// PUBLIC
//
////
public:

  typedef std::vector<double>     NJetScalesType;

  typedef MetaLDA::LDAValuesType  LDAValuesType;

  typedef MetaLDA::LDAMatrixType  LDAMatrixType;

  ////
  //
  // Constructors & Destructor
  //
  ////
  MetaNJetLDA( void );

  MetaNJetLDA( const char *_headerName );

  MetaNJetLDA( const MetaNJetLDA & _metaNJetLDA );

  MetaNJetLDA( const NJetScalesType & _zeroScales,
      const NJetScalesType & _firstScales,
      const NJetScalesType & _secondScales,
      const NJetScalesType & _ridgeScales,
      const LDAValuesType & _ldaValues,
      const LDAMatrixType & _ldaMatrix );

  ~MetaNJetLDA( void );

  void  PrintInfo( void ) const;

  void  CopyInfo( const MetaNJetLDA & _lda );

  void  Clear( void );

  bool  InitializeEssential(
    const NJetScalesType & _zeroScales,
    const NJetScalesType & _firstScales,
    const NJetScalesType & _secondScales,
    const NJetScalesType & _ridgeScales,
    const LDAValuesType & _ldaValues,
    const LDAMatrixType & _ldaMatrix );

  //
  //
  //
  void  SetZeroScales( const NJetScalesType & _zeroScales );
  const NJetScalesType & GetZeroScales( void ) const;

  void  SetFirstScales( const NJetScalesType & _firstScales );
  const NJetScalesType & GetFirstScales( void ) const;

  void  SetSecondScales( const NJetScalesType & _secondScales );
  const NJetScalesType & GetSecondScales( void ) const;

  void  SetRidgeScales( const NJetScalesType & _ridgeScales );
  const NJetScalesType & GetRidgeScales( void ) const;

  //
  //
  //
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

  NJetScalesType  m_ZeroScales;
  NJetScalesType  m_FirstScales;
  NJetScalesType  m_SecondScales;
  NJetScalesType  m_RidgeScales;

  LDAValuesType   m_ZeroScalesTmp;
  LDAValuesType   m_FirstScalesTmp;
  LDAValuesType   m_SecondScalesTmp;
  LDAValuesType   m_RidgeScalesTmp;

  void  M_Destroy( void );

  void  M_SetupReadFields( void );

  void  M_SetupWriteFields( void );

  bool  M_Read( void );

};

}

}

#endif
