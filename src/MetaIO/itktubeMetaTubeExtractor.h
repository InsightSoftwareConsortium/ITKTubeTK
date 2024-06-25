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

#ifndef __itktubeMetaTubeExtractor_h
#define __itktubeMetaTubeExtractor_h

#include <metaTypes.h>
#include <metaForm.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#ifndef METAIO_STREAM
#define METAIO_STREAM std
#endif

namespace itk
{

namespace tube
{

/**
*
* Reads and Writes MetaTubeExtractor Files, typically designated .mtp files
*
* \author Stephen R. Aylward
*
* \date August 29, 1999
*
*/
class MetaTubeExtractor : public MetaForm
{
public:

  typedef std::vector< double >   ValueListType;

  typedef vnl_vector< double >    VectorType;

  typedef vnl_matrix< double >    MatrixType;

  MetaTubeExtractor( void );

  MetaTubeExtractor( const char * _headerName );

  MetaTubeExtractor( const MetaTubeExtractor & _metaTubeExtractor );

  ~MetaTubeExtractor( void );

  virtual void  PrintInfo( void ) const;

  void SetGeneralProperties( double _dataMin, double _dataMax,
    const VectorType & _tubeColor );

  void SetRidgeProperties( double _ridgeScale,
    double _ridgeScaleKernelExtent,
    bool _ridgeDynamicScale,
    bool _ridgeDynamicStepSize,
    double _ridgeStepX,
    double _ridgeMaxTangentChange,
    double _ridgeMaxXChange,
    double _ridgeMinRidgeness,
    double _ridgeMinRidgenessStart,
    double _ridgeMinRoundness,
    double _ridgeMinRoundnessStart,
    double _ridgeMinCurvature,
    double _ridgeMinCurvatureStart,
    double _ridgeMinLevelness,
    double _ridgeMinLevelnessStart,
    int _ridgeMaxRecoveryAttempts );

  void SetRadiusProperties( double _radiusStart,
    double _radiusMin,
    double _radiusMax,
    double _radiusMinMedialness,
    double _radiusMinMedialnessStart );

  double GetDataMin( void ) const;
  double GetDataMax( void ) const;

  VectorType GetTubeColor( void ) const;

  double GetRidgeScale( void ) const;
  double GetRidgeScaleKernelExtent( void ) const;
  bool GetRidgeDynamicScale( void ) const;
  bool GetRidgeDynamicStepSize( void ) const;
  double GetRidgeStepX( void ) const;
  double GetRidgeMaxTangentChange( void ) const;
  double GetRidgeMaxXChange( void ) const;
  double GetRidgeMinRidgeness( void ) const;
  double GetRidgeMinRidgenessStart( void ) const;
  double GetRidgeMinRoundness( void ) const;
  double GetRidgeMinRoundnessStart( void ) const;
  double GetRidgeMinCurvature( void ) const;
  double GetRidgeMinCurvatureStart( void ) const;
  double GetRidgeMinLevelness( void ) const;
  double GetRidgeMinLevelnessStart( void ) const;
  int GetRidgeMaxRecoveryAttempts( void ) const;

  double GetRadiusStart( void ) const;
  double GetRadiusMin( void ) const;
  double GetRadiusMax( void ) const;
  double GetRadiusMinMedialness( void ) const;
  double GetRadiusMinMedialnessStart( void ) const;

  using MetaForm::CopyInfo;
  virtual void  CopyInfo( const MetaTubeExtractor & _tubeExtractor );

  virtual void  Clear( void );

  bool  InitializeEssential( void );

  virtual bool CanRead( const char * _headerName = NULL ) const;
  virtual bool Read( const char * _headerName = NULL );
  virtual bool CanReadStream( METAIO_STREAM::ifstream * _stream ) const;

  virtual bool ReadStream( METAIO_STREAM::ifstream * _stream );

  virtual bool Write( const char *_headName = NULL );

  virtual bool WriteStream( METAIO_STREAM::ofstream * _stream );

protected:

  void  M_Destroy( void );

  void  M_SetupReadFields( void );

  void  M_SetupWriteFields( void );

  bool  M_Read( void );

  double       m_DataMin;
  double       m_DataMax;
  VectorType   m_TubeColor;

  double       m_RidgeScale;
  double       m_RidgeScaleKernelExtent;
  bool         m_RidgeDynamicScale;
  bool         m_RidgeDynamicStepSize;
  double       m_RidgeStepX;
  double       m_RidgeMaxTangentChange;
  double       m_RidgeMaxXChange;
  double       m_RidgeMinRidgeness;
  double       m_RidgeMinRidgenessStart;
  double       m_RidgeMinRoundness;
  double       m_RidgeMinRoundnessStart;
  double       m_RidgeMinCurvature;
  double       m_RidgeMinCurvatureStart;
  double       m_RidgeMinLevelness;
  double       m_RidgeMinLevelnessStart;
  int          m_RidgeMaxRecoveryAttempts;

  double       m_RadiusStart;
  double       m_RadiusMin;
  double       m_RadiusMax;
  double       m_RadiusMinMedialness;
  double       m_RadiusMinMedialnessStart;


}; // End class MetaTubeExtractor

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeMetaTubeExtractor_h )
