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
#ifndef __itkTubeLDAGenerator_h
#define __itkTubeLDAGenerator_h

#include <vector>

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include "itkImage.h"

namespace itk
{

namespace tube
{

template< class ImageT, class LabelmapT >
class LDAGenerator : public Object
{
public:

  typedef LDAGenerator                         Self;
  typedef Object                               Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkTypeMacro( LDAGenerator, Object );

  itkNewMacro( Self );

  //
  // Custom Typedefs
  //
  typedef ImageT                                        ImageType;
  typedef std::vector< typename ImageType::Pointer >    ImageListType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    ImageT::ImageDimension );

  typedef LabelmapT                            MaskImageType;

  typedef double                               FeatureType;
  typedef vnl_vector< FeatureType >            FeatureVectorType;

  typedef int                                  ObjectIdType;
  typedef std::vector< ObjectIdType >          ObjectIdListType;

  typedef vnl_vector< double >                 ObjectMeanType;
  typedef std::vector< ObjectMeanType >        ObjectMeanListType;

  typedef vnl_matrix< double >                 ObjectCovarianceType;
  typedef std::vector< ObjectCovarianceType >  ObjectCovarianceListType;

  typedef vnl_vector< double >                 LDAValuesType;
  typedef vnl_vector< double >                 LDAVectorType;
  typedef vnl_matrix< double >                 LDAMatrixType;

  typedef itk::Image< float, ImageDimension >  LDAImageType;

  //
  // Methods
  //
  void SetFeatureImage( typename ImageType::Pointer img );
  void AddFeatureImage( typename ImageType::Pointer img );
  virtual typename ImageType::Pointer GetFeatureImage( unsigned int num );
  virtual unsigned int GetNumberOfFeatureImages( void );

  virtual unsigned int GetNumberOfFeatures( void );

  void             SetObjectId( ObjectIdType objectId );
  void             AddObjectId( ObjectIdType objectId );
  ObjectIdType     GetObjectId( int num = 0 );

  unsigned int     GetNumberOfObjects( void );

  ObjectMeanType       * GetObjectMean( ObjectIdType objectId );
  ObjectCovarianceType * GetObjectCovariance( ObjectIdType objectId );

  ObjectMeanType       * GetGlobalMean( void );
  ObjectCovarianceType * GetGlobalCovariance( void );

  itkSetObjectMacro( Labelmap, MaskImageType );
  itkGetObjectMacro( Labelmap, MaskImageType );

  unsigned int    GetNumberOfLDA( void );

  LDAVectorType   GetLDAVector( unsigned int ldaNum );
  double          GetLDAValue( unsigned int ldaNum );
  LDAMatrixType & GetLDAMatrix( void );
  LDAValuesType & GetLDAValues( void );
  void            SetLDAVector( unsigned int ldaNum,
                    const LDAVectorType & vec );
  void            SetLDAValue( unsigned int ldaNum, double value );
  void            SetLDAMatrix( const LDAMatrixType & mat );
  void            SetLDAValues( const LDAValuesType & values );

  typename LDAImageType::Pointer GetLDAImage( unsigned int ldaNum );

  itkSetMacro( PerformLDA, bool );
  itkGetMacro( PerformLDA, bool );
  itkSetMacro( PerformPCA, bool );
  itkGetMacro( PerformPCA, bool );

  void SetProgressProcessInformation( void * processInfo, double fraction,
    double start );

  virtual void Update( void );
  virtual void UpdateLDAImages( void );

protected:

  LDAGenerator( void );
  virtual ~LDAGenerator( void );

  typedef ContinuousIndex< double, ImageDimension > ContinuousIndexType;

  virtual LDAValuesType GetFeatureVector( const ContinuousIndexType &
    indx );

  virtual void GenerateStatistics( void );
  virtual void GenerateLDA( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

  FeatureVectorType               m_FeatureVector;

private:

  LDAGenerator( const Self & );          // Purposely not implemented
  void operator = ( const Self & );      // Purposely not implemented

  bool                            m_PerformLDA;
  bool                            m_PerformPCA;

  //  Data
  ImageListType                   m_FeatureImageList;

  typename MaskImageType::Pointer m_Labelmap;

  ObjectIdListType                m_ObjectIdList;
  ObjectMeanListType              m_ObjectMeanList;
  ObjectCovarianceListType        m_ObjectCovarianceList;

  ObjectMeanType                  m_GlobalMean;
  ObjectCovarianceType            m_GlobalCovariance;

  unsigned int                    m_NumberOfLDA;
  LDAMatrixType                   m_LDAMatrix;
  LDAValuesType                   m_LDAValues;

  typename LDAImageType::Pointer  m_LDAImage;

  void                          * m_ProgressProcessInfo;
  double                          m_ProgressFraction;
  double                          m_ProgressStart;


};

}

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeLDAGenerator.txx"
#endif

#endif
