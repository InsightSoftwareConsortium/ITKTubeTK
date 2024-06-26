/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeMetaRidgeSeed_h
#define __itktubeMetaRidgeSeed_h


#include "itktubeMetaLDA.h"

#include <metaTypes.h>
#include <metaForm.h>

#ifndef METAIO_STREAM
#define METAIO_STREAM std
#endif

namespace itk
{

namespace tube
{

/**
*
* Reads and Writes MetaRidgeSeed Files, typically designated .mnda files
*
* REQUIRED: RidgeSeedGenerator instance
*
* \author Stephen R. Aylward
*
* \date August 29, 1999
*
*/
class MetaRidgeSeed : public MetaLDA
{
public:

  typedef std::vector< double >   RidgeSeedScalesType;

  typedef MetaLDA::LDAValuesType  LDAValuesType;

  typedef MetaLDA::LDAMatrixType  LDAMatrixType;

  MetaRidgeSeed( void );

  MetaRidgeSeed( const char * _headerName );

  MetaRidgeSeed( const MetaRidgeSeed & _metaRidgeSeed );

  MetaRidgeSeed(
    const RidgeSeedScalesType & _ridgeScales,
    const bool _useIntensityOnly,
    const bool _useFeatureMath,
    const LDAValuesType & _ldaValues,
    const LDAMatrixType & _ldaMatrix,
    const ValueListType & _inputWhitenMeans,
    const ValueListType & _inputWhitenStdDevs,
    const ValueListType & _outputWhitenMeans,
    const ValueListType & _outputWhitenStdDevs,
    const std::string & _pdfFileName );

  ~MetaRidgeSeed( void );

  virtual void PrintInfo( void ) const;

  using MetaLDA::CopyInfo;
  virtual void CopyInfo( const MetaRidgeSeed & _lda );

  virtual void Clear( void );

  bool InitializeEssential(
    const RidgeSeedScalesType & _ridgeScales,
    const bool _useIntensityOnly,
    const bool _useFeatureMath,
    const LDAValuesType & _ldaValues,
    const LDAMatrixType & _ldaMatrix,
    const ValueListType & _inputWhitenMeans,
    const ValueListType & _inputWhitenStdDevs,
    const ValueListType & _outputWhitenMeans,
    const ValueListType & _outputWhitenStdDevs,
    const std::string & _pdfFileName );

  void  SetRidgeSeedScales( const RidgeSeedScalesType & _ridgeScales );
  const RidgeSeedScalesType & GetRidgeSeedScales( void ) const;

  void SetUseIntensityOnly( bool _useIntensityOnly );
  bool GetUseIntensityOnly( void ) const;

  void SetUseFeatureMath( bool _useFeatureMath );
  bool GetUseFeatureMath( void ) const;

  void  SetPDFFileName( const std::string & _pdfFileName );
  const std::string & GetPDFFileName( void ) const;

  void SetRidgeId( int _ridgeId );
  int  GetRidgeId( void ) const;

  void SetBackgroundId( int _backgroundId );
  int  GetBackgroundId( void ) const;

  void SetUnknownId( int _unknownId );
  int  GetUnknownId( void ) const;

  void   SetSeedTolerance( double _seedTolerance );
  double GetSeedTolerance( void ) const;

  void SetSkeletonize( bool _skeletonize );
  bool GetSkeletonize( void ) const;

  virtual bool CanRead( const char * _headerName = NULL ) const;

  virtual bool Read( const char * _headerName = NULL );

  virtual bool CanReadStream( METAIO_STREAM::ifstream * _stream ) const;

  virtual bool ReadStream( METAIO_STREAM::ifstream * _stream );

  virtual bool Write( const char * _headName = NULL );

  virtual bool WriteStream( METAIO_STREAM::ofstream * _stream );

protected:
  void  M_Destroy( void );

  void  M_SetupReadFields( void );

  void  M_SetupWriteFields( void );

  bool  M_Read( void );

  bool   m_UseIntensityOnly;
  bool   m_UseFeatureMath;

  int    m_RidgeId;
  int    m_BackgroundId;
  int    m_UnknownId;
  double m_SeedTolerance;
  bool   m_Skeletonize;

  RidgeSeedScalesType  m_RidgeSeedScales;

  std::string          m_PDFFileName;

}; // End class MetaRidgeSeed

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeMetaRidgeSeed_h )
