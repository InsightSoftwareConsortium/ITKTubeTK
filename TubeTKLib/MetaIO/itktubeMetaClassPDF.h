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

#ifndef __itktubeMetaClassPDF_h
#define __itktubeMetaClassPDF_h

#include "metaImage.h"

#include <metaForm.h>

namespace itk
{

namespace tube
{

/**
*
* Reads and Writes MetaClassPDF Files, typically designated .mnda files
*
* \author Stephen R. Aylward
*
* \date August 29, 2013
*
*/
class MetaClassPDF : private MetaImage
{
public:

  typedef std::vector< int >                    VectorIntType;
  typedef std::vector< unsigned int >           VectorUIntType;
  typedef std::vector< double >                 VectorDoubleType;

  MetaClassPDF( void );

  MetaClassPDF( const char * _headerName );

  MetaClassPDF( const MetaClassPDF & _metaPDF );

  MetaClassPDF( unsigned int _nFeatures,
    const VectorUIntType & _nBinsPerFeature,
    const VectorDoubleType & _binMin,
    const VectorDoubleType & _binSize,
    float * _elementData = NULL );

  MetaClassPDF( unsigned int _x,
    unsigned int _y,
    double _binMinX,
    double _binMinY,
    double _binSizeX,
    double _binSizeY,
    float * _elementData = NULL );

  MetaClassPDF( unsigned int _x,
    unsigned int _y,
    unsigned int _z,
    double _binMinX,
    double _binMinY,
    double _binMinZ,
    double _binSizeX,
    double _binSizeY,
    double _binSizeZ,
    float * _elementData = NULL );

  ~MetaClassPDF( void );

  virtual void PrintInfo( void ) const;

  virtual void CopyInfo( const MetaObject * _pdf );

  virtual void Clear( void );

  virtual bool InitializeEssential( unsigned int _nFeatures,
    const VectorUIntType & _nBinsPerFeature,
    const VectorDoubleType & _binMin,
    const VectorDoubleType & _binSize,
    float * _elementData = NULL );

  void         SetNumberOfFeatures( unsigned int _nFeatures );
  unsigned int GetNumberOfFeatures( void ) const;

  void         SetNumberOfBinsPerFeature( const VectorUIntType & _nBins );
  const VectorUIntType & GetNumberOfBinsPerFeature( void ) const;

  void          SetBinMin( const VectorDoubleType & _binMin );
  const VectorDoubleType & GetBinMin( void ) const;

  void          SetBinSize( const VectorDoubleType & _binSize );
  const VectorDoubleType & GetBinSize( void ) const;

  void          SetPDF( float * _pdfData );
  float *       GetPDF( void );  // Data is freed when reader is destroyed
  float *       ExportPDF( void );  // Data persists when reader destroyed

  void          SetObjectId( const VectorIntType & _objectIds );
  const VectorIntType & GetObjectId( void ) const;

  void          SetObjectPDFWeight( const VectorDoubleType &
                  _objectWeights );
  const VectorDoubleType & GetObjectPDFWeight( void ) const;

  void          SetVoidId( int _voidId );
  int           GetVoidId( void ) const;

  void          SetErodeRadius( unsigned int _ErodeRadius );
  unsigned int  GetErodeRadius( void ) const;
  void          SetHoleFillIterations( unsigned int _HoleFillIterations );
  unsigned int  GetHoleFillIterations( void ) const;
  void          SetProbabilityImageSmoothingStandardDeviation(
                  double _ProbabilityImageSmoothingStandardDeviation );
  double         GetProbabilityImageSmoothingStandardDeviation( void ) const;
  void          SetHistogramSmoothingStandardDeviation(
                  double _HistogramSmoothingStandardDeviation );
  double         GetHistogramSmoothingStandardDeviation( void ) const;
  void          SetOutlierRejectPortion( double _OutlierRejectPortion );
  double         GetOutlierRejectPortion( void ) const;
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

  virtual bool CanReadStream( std::ifstream * _stream ) const;

  virtual bool ReadStream( std::ifstream * _stream );

  virtual bool Write( const char * _headerName = NULL );

  virtual bool WriteStream( std::ofstream * _stream );

protected:

  virtual void M_SetupReadFields( void );

  virtual void M_SetupWriteFields( void );

  virtual bool M_Read( void );

  mutable VectorUIntType   m_NumberOfBinsPerFeature;
  mutable VectorDoubleType m_BinMin;
  mutable VectorDoubleType m_BinSize;

  VectorIntType        m_ObjectId;
  VectorDoubleType     m_ObjectPDFWeight;
  int                  m_VoidId;
  unsigned int         m_ErodeRadius;
  unsigned int         m_HoleFillIterations;
  double               m_ProbabilityImageSmoothingStandardDeviation;
  double               m_HistogramSmoothingStandardDeviation;
  double               m_OutlierRejectPortion;
  bool                 m_Draft;
  bool                 m_ReclassifyObjectLabels;
  bool                 m_ReclassifyNotObjectLabels;
  bool                 m_ForceClassification;

}; // End class MetaRidgeSeed

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeMetaClassPDF_h )
