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

#ifndef __itkTubeMetaLDA_h
#define __itkTubeMetaLDA_h

#include <metaForm.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

namespace itk
{

namespace tube
{

/**
 * \brief    Reads and writes MetaLDA files, typically designated .mlda files.
 * \pre      itkTubeLDAGenerator instance.
 *
 * \author   Stephen R. Aylward
 * \date     August 29, 1999
 *
 * \ingroup  IO
 */
class METAIO_EXPORT MetaLDA : public MetaForm
{
public:

  typedef std::vector< double >  ValueListType;
  typedef vnl_vector< double >   LDAValuesType;
  typedef vnl_matrix< double >   LDAMatrixType;

  MetaLDA( void );

  MetaLDA( const char * headerName );

  MetaLDA( const MetaLDA & metaLDA );

  MetaLDA( const LDAValuesType & ldaValues,
           const LDAMatrixType & ldaMatrix,
           const ValueListType & whitenMeans,
           const ValueListType & whitenStdDevs );

  ~MetaLDA( void );

  virtual void PrintInfo( void ) const;

  using MetaForm::CopyInfo;
  virtual void CopyInfo( const MetaLDA & lda );

  virtual void Clear( void );

  bool InitializeEssential( const LDAValuesType & ldaValues,
                            const LDAMatrixType & ldaMatrix,
                            const ValueListType & whitenMeans,
                            const ValueListType & whitenStdDevs );

  void SetLDAValues( const LDAValuesType & ldaValues );

  const LDAValuesType & GetLDAValues( void ) const;

  void SetLDAMatrix( const LDAMatrixType & ldaMatrix );

  const LDAMatrixType & GetLDAMatrix( void ) const;

  void SetWhitenMeans( const ValueListType & whitenMeans );

  const ValueListType & GetWhitenMeans( void ) const;

  void SetWhitenStdDevs( const ValueListType & whitenStdDevs );

  const ValueListType & GetWhitenStdDevs( void ) const;

  virtual bool CanRead( const char * headerName = NULL ) const;

  virtual bool Read( const char * headerName = NULL );

  virtual bool CanReadStream( METAIO_STREAM::ifstream * stream ) const;

  virtual bool ReadStream( METAIO_STREAM::ifstream * stream );

  virtual bool Write( const char * headerName = NULL );

  virtual bool WriteStream( METAIO_STREAM::ofstream * stream );

protected:

  void M_Destroy( void );

  void M_SetupReadFields( void );

  void M_SetupWriteFields( void );

  bool M_Read( void );

  LDAValuesType  m_LDAValues;

  ValueListType  m_WhitenMeans;
  ValueListType  m_WhitenStdDevs;

  LDAMatrixType  m_LDAMatrix;

}; // End class MetaLDA

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itkTubeMetaLDA_h)
