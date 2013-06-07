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

#ifndef __itkTubeMetaTubeParams_h
#define __itkTubeMetaTubeParams_h

#include <metaForm.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

namespace itk
{

namespace tube
{

/**
*
* Reads and Writes MetaTubeParams Files, typically designated .mtp files
*
* \author Stephen R. Aylward
*
* \date August 29, 1999
*
*/
class METAIO_EXPORT MetaTubeParams
: public MetaForm
{
public:

  typedef std::vector< double >   ValueListType;

  typedef vnl_vector< double >    VectorType;

  typedef vnl_matrix< double >    MatrixType;

  MetaTubeParams( void );

  MetaTubeParams( const char * _headerName );

  MetaTubeParams( const MetaTubeParams & _metaTubeParams );

  ~MetaTubeParams( void );

  virtual void  PrintInfo( void ) const;

  void SetSeedParams( const VectorType & _seedScales,
    double _seedIntensityMin, double _seedIntensityMax,
    double _seedPercentile );

  void SetTubeParams( double _tubeIntensityMin, double _tubeIntensityMax,
    bool _tubeBright, const VectorType & _tubeColor );

  void SetTubeRidgeParams( double _tubeRidgeScale,
    double _tubeRidgeScaleExtent,
    bool _tubeRidgeDynamicScale, double _tubeRidgeStepX,
    double _tubeRidgeThresholdTangentChange,
    double _tubeRidgeThresholdXChange,
    double _tubeRidgeThresholdRidgeness,
    double _tubeRidgeThresholdRidgenessStart,
    double _tubeRidgeThresholdRoundness,
    double _tubeRidgeThresholdRoundnessStart,
    double _tubeRidgeCurvatureMax,
    double _tubeRidgeThresholdCurvature,
    double _tubeRidgeThresholdCurvatureStart,
    double _tubeRidgeThresholdLinearity,
    double _tubeRidgeThresholdLinearityStart,
    int _tubeRidgeRecoveryMax );

  void SetTubeRadiusParams( double _tubeRadiusStart,
    double _tubeRadiusMin, double _tubeRadiusMax,
    double _tubeRadiusThresholdMedialness,
    double _tubeRadiusThresholdMedialnessStart );

  VectorType GetSeedScales( void ) const;
  double GetSeedIntensityMin( void ) const;
  double GetSeedIntensityMax( void ) const;
  double GetSeedIntensityPercentile( void ) const;

  double GetTubeIntensityMin( void ) const;
  double GetTubeIntensityMax( void ) const;
  bool GetTubeBright( void ) const;
  VectorType GetTubeColor( void ) const;

  double GetTubeRidgeScale( void ) const;
  double GetTubeRidgeScaleExtent( void ) const;
  bool GetTubeRidgeDynamicScale( void ) const;
  double GetTubeRidgeStepX( void ) const;
  double GetTubeRidgeThresholdTangentChange( void ) const;
  double GetTubeRidgeThresholdXChange( void ) const;
  double GetTubeRidgeThresholdRidgeness( void ) const;
  double GetTubeRidgeThresholdRidgenessStart( void ) const;
  double GetTubeRidgeThresholdRoundness( void ) const;
  double GetTubeRidgeThresholdRoundnessStart( void ) const;
  double GetTubeRidgeCurvatureMax( void ) const;
  double GetTubeRidgeThresholdCurvature( void ) const;
  double GetTubeRidgeThresholdCurvatureStart( void ) const;
  double GetTubeRidgeThresholdLinearity( void ) const;
  double GetTubeRidgeThresholdLinearityStart( void ) const;
  int GetTubeRidgeRecoveryMax( void ) const;

  double GetTubeRadiusStart( void ) const;
  double GetTubeRadiusMin( void ) const;
  double GetTubeRadiusMax( void ) const;
  double GetTubeRadiusThresholdMedialness( void ) const;
  double GetTubeRadiusThresholdMedialnessStart( void ) const;

  void  CopyInfo( const MetaTubeParams & _tubeParams );

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

  VectorType   m_SeedScales;
  double       m_SeedIntensityMin;
  double       m_SeedIntensityMax;
  double       m_SeedIntensityPercentile;

  double       m_TubeIntensityMin;
  double       m_TubeIntensityMax;
  bool         m_TubeBright;
  VectorType   m_TubeColor;

  double       m_TubeRidgeScale;
  double       m_TubeRidgeScaleExtent;
  bool         m_TubeRidgeDynamicScale;
  double       m_TubeRidgeStepX;
  double       m_TubeRidgeThresholdTangentChange;
  double       m_TubeRidgeThresholdXChange;
  double       m_TubeRidgeThresholdRidgeness;
  double       m_TubeRidgeThresholdRidgenessStart;
  double       m_TubeRidgeThresholdRoundness;
  double       m_TubeRidgeThresholdRoundnessStart;
  double       m_TubeRidgeCurvatureMax;
  double       m_TubeRidgeThresholdCurvature;
  double       m_TubeRidgeThresholdCurvatureStart;
  double       m_TubeRidgeThresholdLinearity;
  double       m_TubeRidgeThresholdLinearityStart;
  int          m_TubeRidgeRecoveryMax;

  double       m_TubeRadiusStart;
  double       m_TubeRadiusMin;
  double       m_TubeRadiusMax;
  double       m_TubeRadiusThresholdMedialness;
  double       m_TubeRadiusThresholdMedialnessStart;


}; // End class MetaTubeParams

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itkTubeMetaTubeParams_h)
