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

#ifndef __itktubeMetaPDF_H
#define __itktubeMetaPDF_H

#include "metaImage.h"

#include <metaForm.h>

namespace itk
{

namespace tube
{

/**
*
* Reads and Writes MetaPDF Files, typically designated .mnda files
*
* \author Stephen R. Aylward
*
* \date August 29, 2013
*
*/
class MetaPDF : private MetaImage
{
public:

  MetaPDF( void );

  MetaPDF( const char * _headerName );

  MetaPDF( const MetaPDF & _metaPDF );

  MetaPDF( int _nFeatures,
    const std::vector< int > & _nBinsPerFeature,
    const std::vector< float > & _binMin,
    const std::vector< float > & _binSize,
    float * _elementData = NULL );

  MetaPDF( int _x, int _y,
    float _binMinX, float _binMinY,
    float _binSizeX, float _binSizeY,
    float * _elementData = NULL );

  MetaPDF( int _x, int _y, int _z,
    float _binMinX, float _binMinY, float _binMinZ,
    float _binSizeX, float _binSizeY, float _binSizeZ,
    float * _elementData = NULL );

  ~MetaPDF( void );

  virtual void PrintInfo( void ) const;

  virtual void CopyInfo( const MetaPDF & _pdf );

  virtual void Clear( void );

  virtual bool InitializeEssential( int _nFeatures,
    const std::vector< int > & _nBinsPerFeature,
    const std::vector< float > & _binMin,
    const std::vector< float > & _binSize,
    float * _elementData = NULL );

  void         SetNumberOfFeatures( int _nFeatures );
  int          GetNumberOfFeatures( void ) const;

  void         SetNumberOfBinsPerFeature( const std::vector< int > &
                 _nBins );
  const std::vector< int > & GetNumberOfBinsPerFeature( void ) const;

  void          SetBinMin( const std::vector< float > & _binMin );
  const std::vector< float > & GetBinMin( void ) const;

  void          SetBinSize( const std::vector< float > & _binSize );
  const std::vector< float > & GetBinSize( void ) const;

  void          SetPDF( float * _pdfData );
  float *       GetPDF( void );

  void          SetObjectId( const std::vector< int > & _objectIds );
  const std::vector< int > & GetObjectId( void ) const;

  void          SetObjectPDFWeight( const std::vector< float > &
                  _objectWeights );
  const std::vector< float > & GetObjectPDFWeight( void ) const;

  void          SetVoidId( int _voidId );
  int           GetVoidId( void ) const;

  void          SetErodeRadius( int _ErodeRadius );
  int           GetErodeRadius( void ) const;
  void          SetHoleFillIterations( int _HoleFillIterations );
  int           GetHoleFillIterations( void ) const;
  void          SetProbabilityImageSmoothingStandardDeviation(
                  float _ProbabilityImageSmoothingStandardDeviation );
  float         GetProbabilityImageSmoothingStandardDeviation( void ) const;
  void          SetHistogramSmoothingStandardDeviation(
                  float _HistogramSmoothingStandardDeviation );
  float         GetHistogramSmoothingStandardDeviation( void ) const;
  void          SetOutlierRejectPortion( float _OutlierRejectPortion );
  float         GetOutlierRejectPortion( void ) const;
  void          SetDraft( bool _Draft );
  bool          GetDraft( void ) const;
  void          SetReclassifyObjectLabels( bool _ReclassifyObjectLabels );
  bool          GetReclassifyObjectLabels( void ) const;
  void          SetReclassifyNotObjectLabels(
                  bool _ReclassifyNotObjectLabels );
  bool          GetReclassifyNotObjectLabels( void ) const;
  void          SetForceClassification( bool _ForceClassification );
  bool          GetForceClassification( void ) const;

  virtual bool CanRead( const char * _headerName = NULL ) const;

  virtual bool Read( const char * _headerName = NULL );

  virtual bool CanReadStream( METAIO_STREAM::ifstream * _stream ) const;

  virtual bool ReadStream( METAIO_STREAM::ifstream * _stream );

  virtual bool Write( const char * _headerName = NULL );

  virtual bool WriteStream( METAIO_STREAM::ofstream * _stream );

protected:

  virtual void M_SetupReadFields( void );

  virtual void M_SetupWriteFields( void );

  virtual bool M_Read( void );

  mutable std::vector< int >   m_NumberOfBinsPerFeature;
  mutable std::vector< float > m_BinMin;
  mutable std::vector< float > m_BinSize;

  std::vector< int >   m_ObjectId;
  std::vector< float > m_ObjectPDFWeight;
  int                  m_VoidId;
  int                  m_ErodeRadius;
  int                  m_HoleFillIterations;
  float                m_ProbabilityImageSmoothingStandardDeviation;
  float                m_HistogramSmoothingStandardDeviation;
  float                m_OutlierRejectPortion;
  bool                 m_Draft;
  bool                 m_ReclassifyObjectLabels;
  bool                 m_ReclassifyNotObjectLabels;
  bool                 m_ForceClassification;

}; // End class MetaRidgeSeed

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeMetaPDF_H)
