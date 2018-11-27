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

#ifndef __itktubeMetaLDA_h
#define __itktubeMetaLDA_h

#include <metaForm.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_c_vector.h>

namespace itk
{

namespace tube
{

/**
 * \brief    Reads and writes MetaLDA files, typically designated .mlda files.
 * \pre      LDAGenerator instance.
 *
 * \author   Stephen R. Aylward
 * \date     August 29, 1999
 *
 * \ingroup  IO
 */
class MetaLDA : public MetaForm
{
public:

  typedef std::vector< double >  ValueListType;
  typedef vnl_vector< double >   LDAValuesType;
  typedef vnl_matrix< double >   LDAMatrixType;

  MetaLDA( void );

  MetaLDA( const char * headerName );

  MetaLDA( const MetaLDA & metaLDA );

  MetaLDA( unsigned int _numberOfPCABasisToUseAsFeatures,
    unsigned int _numberOfLDABasisToUseAsFeatures,
    const LDAValuesType & _ldaValues,
    const LDAMatrixType & _ldaMatrix,
    const ValueListType & _inputWhitenMeans,
    const ValueListType & _inputWhitenStdDevs,
    const ValueListType & _outputWhitenMeans,
    const ValueListType & _outputWhitenStdDevs );

  ~MetaLDA( void );

  virtual void PrintInfo( void ) const;

  using MetaForm::CopyInfo;
  virtual void CopyInfo( const MetaLDA & lda );

  virtual void Clear( void );

  bool  InitializeEssential( unsigned int _numberOfPCABasisToUseAsFeatures,
    unsigned int _numberOfLDABasisToUseAsFeatures,
    const LDAValuesType & _ldaValues,
    const LDAMatrixType & _ldaMatrix,
    const ValueListType & _inputWhitenMeans,
    const ValueListType & _inputWhitenStdDevs,
    const ValueListType & _outputWhitenMeans,
    const ValueListType & _outputWhitenStdDevs );

  void SetNumberOfPCABasisToUseAsFeatures( unsigned int
    numberOfPCABasisToUseAsFeatures );
  unsigned int GetNumberOfPCABasisToUseAsFeatures( void ) const;

  void SetNumberOfLDABasisToUseAsFeatures( unsigned int
    numberOfLDABasisToUseAsFeatures );
  unsigned int GetNumberOfLDABasisToUseAsFeatures( void ) const;

  void SetLDAValues( const LDAValuesType & ldaValues );
  const LDAValuesType & GetLDAValues( void ) const;

  void SetLDAMatrix( const LDAMatrixType & ldaMatrix );
  const LDAMatrixType & GetLDAMatrix( void ) const;

  void SetInputWhitenMeans( const ValueListType & whitenMeans );
  const ValueListType & GetInputWhitenMeans( void ) const;

  void SetInputWhitenStdDevs( const ValueListType & whitenStdDevs );
  const ValueListType & GetInputWhitenStdDevs( void ) const;

  void SetOutputWhitenMeans( const ValueListType & whitenMeans );
  const ValueListType & GetOutputWhitenMeans( void ) const;

  void SetOutputWhitenStdDevs( const ValueListType & whitenStdDevs );
  const ValueListType & GetOutputWhitenStdDevs( void ) const;

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

  unsigned int   m_NumberOfPCABasisToUseAsFeatures;

  unsigned int   m_NumberOfLDABasisToUseAsFeatures;

  LDAValuesType  m_LDAValues;

  ValueListType  m_InputWhitenMeans;
  ValueListType  m_InputWhitenStdDevs;

  ValueListType  m_OutputWhitenMeans;
  ValueListType  m_OutputWhitenStdDevs;

  LDAMatrixType  m_LDAMatrix;

}; // End class MetaLDA

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeMetaLDA_h )
