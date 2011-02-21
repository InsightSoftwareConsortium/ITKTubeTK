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
#ifndef __itkTubeLDAGenerator_txx
#define __itkTubeLDAGenerator_txx

#include <limits>

#include "itkTubeLDAGenerator.h"

#include "itkTimeProbesCollectorBase.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "tubeMatrixMath.h"

namespace itk
{

namespace tube
{

template< class ImageT, class LabelmapT >
LDAGenerator< ImageT, LabelmapT >
::LDAGenerator()
{
  m_PerformLDA = true;
  m_PerformPCA = false;

  m_FeatureImageList.clear();

  m_Labelmap = NULL;

  m_ObjectIdList.clear();
  m_ObjectIdList.push_back( 1 );
  m_ObjectMeanList.clear();
  m_ObjectCovarianceList.clear();

  m_LDAValues.set_size( 0 );
  m_LDAMatrix.set_size( 0, 0 );

  m_LDAImage = NULL;

  m_ProgressProcessInfo = NULL;
  m_ProgressFraction = 1.0;
  m_ProgressStart = 0;
}

template< class ImageT, class LabelmapT >
LDAGenerator< ImageT, LabelmapT >
::~LDAGenerator()
{
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::SetFeatureImage( typename ImageType::Pointer img )
{
  m_FeatureImageList.clear();
  m_FeatureImageList.push_back( img );
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::AddFeatureImage( typename ImageType::Pointer img )
{
  m_FeatureImageList.push_back( img );
}

template < class ImageT, class LabelmapT >
typename ImageT::Pointer
LDAGenerator< ImageT, LabelmapT >
::GetFeatureImage( unsigned int num )
{
  if( num < m_FeatureImageList.size() )
    {
    return m_FeatureImageList[num];
    }
  else
    {
    return NULL;
    }
}

template < class ImageT, class LabelmapT >
unsigned int
LDAGenerator< ImageT, LabelmapT >
::GetNumberOfFeatureImages( void )
{
  return m_FeatureImageList.size();
}

template < class ImageT, class LabelmapT >
unsigned int
LDAGenerator< ImageT, LabelmapT >
::GetNumberOfFeatures( void )
{
  return m_FeatureImageList.size();
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::SetObjectId( ObjectIdType objectId )
{
  m_ObjectIdList.clear();
  m_ObjectIdList.push_back( objectId );
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::AddObjectId( ObjectIdType objectId )
{
  m_ObjectIdList.push_back( objectId );
}

template < class ImageT, class LabelmapT >
unsigned int
LDAGenerator< ImageT, LabelmapT >
::GetNumberOfObjects( void )
{
   return m_ObjectIdList.size();
}

template < class ImageT, class LabelmapT >
typename LDAGenerator< ImageT, LabelmapT >::ObjectIdType
LDAGenerator< ImageT, LabelmapT >
::GetObjectId( int num )
{
  if( num < m_ObjectIdList.size() )
    {
    return m_ObjectIdList[ num ];
    }
  else
    {
    // voidId
    return m_ObjectIdList[ m_ObjectIdList.size()-1 ];
    }
}

template < class ImageT, class LabelmapT >
typename LDAGenerator< ImageT, LabelmapT >::ObjectMeanType *
LDAGenerator< ImageT, LabelmapT >
::GetObjectMean( int num )
{
  if( num < m_ObjectIdList.size() )
    {
    return & m_ObjectMeanList[ num ];
    }
  else
    {
    // voidId
    return & m_ObjectMeanList[ m_ObjectIdList.size()-1 ];
    }
}

template < class ImageT, class LabelmapT >
typename LDAGenerator< ImageT, LabelmapT >::ObjectCovarianceType *
LDAGenerator< ImageT, LabelmapT >
::GetObjectCovariance( int num )
{
  if( num < m_ObjectIdList.size() )
    {
    return & m_ObjectCovarianceList[ num ];
    }
  else
    {
    // voidId
    return & m_ObjectCovarianceList[ m_ObjectIdList.size()-1 ];
    }
}

template < class ImageT, class LabelmapT >
unsigned int
LDAGenerator< ImageT, LabelmapT >
::GetNumberOfLDA( void )
{
  return m_LDAValues.size();
}

template < class ImageT, class LabelmapT >
typename LDAGenerator< ImageT, LabelmapT >::LDAVectorType
LDAGenerator< ImageT, LabelmapT >
::GetLDAVector( unsigned int ldaNum )
{
  if( ldaNum < m_LDAValues.size() )
    {
    return m_LDAMatrix.get_column( ldaNum );
    }
  else
    {
    m_FeatureVector.set_size( m_LDAValues.size() );
    m_FeatureVector.fill( 0 );
    return m_FeatureVector;
    }
}

template < class ImageT, class LabelmapT >
double
LDAGenerator< ImageT, LabelmapT >
::GetLDAValue( unsigned int ldaNum )
{
  if( ldaNum < m_LDAValues.size() )
    {
    return m_LDAValues[ ldaNum ];
    }
  else
    {
    return 0;
    }
}

template < class ImageT, class LabelmapT >
typename LDAGenerator< ImageT, LabelmapT >::LDAMatrixType &
LDAGenerator< ImageT, LabelmapT >
::GetLDAMatrix( void )
{
  return m_LDAMatrix;
}

template < class ImageT, class LabelmapT >
typename LDAGenerator< ImageT, LabelmapT >::LDAValuesType &
LDAGenerator< ImageT, LabelmapT >
::GetLDAValues( void )
{
  return m_LDAValues;
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::SetLDAValue( unsigned int ldaNum, double value )
{
  if( ldaNum < m_LDAValues.size() )
    {
    m_LDAValues[ ldaNum ] = value;
    }
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::SetLDAMatrix( const LDAMatrixType & mat )
{
  m_LDAMatrix = mat;
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::SetLDAVector( unsigned int ldaNum, const LDAVectorType & vec )
{
  if( ldaNum < m_LDAValues.size() )
    {
    m_LDAMatrix.set_column( ldaNum, vec );
    }
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::SetLDAValues( const LDAValuesType & values )
{
  m_LDAValues = values;
}

template < class ImageT, class LabelmapT >
typename LDAGenerator< ImageT, LabelmapT >::LDAImageType::Pointer
LDAGenerator< ImageT, LabelmapT >
::GetLDAImage( unsigned int ldaNum )
{
  if( ldaNum < m_LDAValues.size()  )
    {
    itk::TimeProbesCollectorBase timeCollector;

    timeCollector.Start( "GenerateLDAImage" );

    typedef itk::ImageRegionIteratorWithIndex< LDAImageType >
      ImageIteratorType;

    unsigned int numFeatures = this->GetNumberOfFeatures();

    typename LDAImageType::RegionType region;
    region = m_FeatureImageList[0]->GetLargestPossibleRegion();

    m_LDAImage = LDAImageType::New();
    m_LDAImage->SetRegions( region );
    m_LDAImage->CopyInformation( m_FeatureImageList[0] );
    m_LDAImage->Allocate();

    ImageIteratorType itLDAIm( m_LDAImage,
      m_LDAImage->GetLargestPossibleRegion() );

    FeatureVectorType v( numFeatures );
    FeatureVectorType vLDA( numFeatures );
    while( !itLDAIm.IsAtEnd() )
      {
      ContinuousIndexType indx = itLDAIm.GetIndex();
      v = this->GetFeatureVector( indx );

      vLDA = v * m_LDAMatrix;

      itLDAIm.Set( vLDA[ ldaNum ] );
      ++itLDAIm;
      }

    timeCollector.Stop( "GenerateLDAImage" );
    timeCollector.Report();

    return m_LDAImage;
    }
  else
    {
    return NULL;
    }
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::SetProgressProcessInformation( void * processInfo, double fraction,
  double start )
{
  m_ProgressProcessInfo = processInfo;
  m_ProgressFraction = fraction;
  m_ProgressStart = start;
}


template < class ImageT, class LabelmapT >
vnl_vector< double >
LDAGenerator< ImageT, LabelmapT >
::GetFeatureVector( const ContinuousIndexType & indx )
{
  unsigned int numFeatureImages = this->GetNumberOfFeatureImages();

  m_FeatureVector.set_size( numFeatureImages );

  unsigned int vCount = 0;
  typename ImageType::IndexType indxI;
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    indxI[i] = static_cast<int>( indx[i] );
    }
  for( unsigned int i=0; i<numFeatureImages; i++ )
    {
    m_FeatureVector[vCount++] = static_cast< FeatureType >(
      m_FeatureImageList[i]->GetPixel( indxI ) );
    }

  return m_FeatureVector;
}


template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::GenerateStatistics( void )
{
  itk::TimeProbesCollectorBase timeCollector;

  timeCollector.Start( "GenerateStatistics" );

  typedef itk::ImageRegionConstIteratorWithIndex< MaskImageType >
    ConstMaskImageIteratorType;
  ConstMaskImageIteratorType itInMask( m_Labelmap,
    m_Labelmap->GetLargestPossibleRegion() );

  unsigned int numClasses = this->GetNumberOfObjects();
  unsigned int numFeatures = this->GetNumberOfFeatures();

  m_ObjectMeanList.resize( numClasses );
  m_ObjectCovarianceList.resize( numClasses );
  std::vector< unsigned int > countList( numClasses );
  ObjectMeanListType sumList( numClasses );
  ObjectCovarianceListType sumOfSquaresList( numClasses );
  for( unsigned int i=0; i<numClasses; i++ )
    {
    m_ObjectMeanList[i].set_size( numFeatures );
    m_ObjectMeanList[i].fill( 0 );
    m_ObjectCovarianceList[i].set_size( numFeatures, numFeatures );
    m_ObjectCovarianceList[i].fill( 0 );
    countList[i] = 0;
    sumList[i].set_size( numFeatures );
    sumList[i].fill( 0 );
    sumOfSquaresList[i].set_size( numFeatures, numFeatures );
    sumOfSquaresList[i].fill( 0 );
    }

  m_GlobalMean.set_size( numFeatures );
  m_GlobalMean.fill( 0 );
  m_GlobalCovariance.set_size( numFeatures, numFeatures );
  m_GlobalCovariance.fill( 0 );
  unsigned int globalCount = 0;
  ObjectMeanType globalSum( numFeatures );
  globalSum.fill( 0 );
  ObjectCovarianceType globalSumOfSquares( numFeatures, numFeatures );
  globalSumOfSquares.fill( 0 );

  itInMask.GoToBegin();
  while( !itInMask.IsAtEnd() )
    {
    ObjectIdType val = static_cast<ObjectIdType>( itInMask.Get() );
    unsigned int valC = 0;
    bool found = false;
    for( unsigned int c=0; c<numClasses; c++ )
      {
      if( val == m_ObjectIdList[c] )
        {
        valC = c;
        found = true;
        break;
        }
      }

    if( found )
      {
      ContinuousIndexType indx = itInMask.GetIndex();
      LDAValuesType v = this->GetFeatureVector( indx );
      ++globalCount;
      for( unsigned int i=0; i<numFeatures; i++ )
        {
        globalSum[i] += v[i];
        for( unsigned int j=0; j<numFeatures; j++ )
          {
          globalSumOfSquares[i][j] += v[i] * v[j];
          }
        }
      }
    ++itInMask;
    }

  for( unsigned int i=0; i<numFeatures; i++ )
    {
    m_GlobalMean[i] = globalSum[i] / globalCount;
    }
  for( unsigned int i=0; i<numFeatures; i++ )
    {
    for( unsigned int j=0; j<numFeatures; j++ )
      {
      m_GlobalCovariance[i][j] = ( globalSumOfSquares[i][j]
        / globalCount ) - ( m_GlobalMean[i] * m_GlobalMean[j] );
      }
    }

  itInMask.GoToBegin();
  while( !itInMask.IsAtEnd() )
    {
    ObjectIdType val = static_cast<ObjectIdType>( itInMask.Get() );
    unsigned int valC = 0;
    bool found = false;
    for( unsigned int c=0; c<numClasses; c++ )
      {
      if( val == m_ObjectIdList[c] )
        {
        valC = c;
        found = true;
        break;
        }
      }

    if( found )
      {
      ContinuousIndexType indx = itInMask.GetIndex();
      LDAValuesType v = this->GetFeatureVector( indx );
      for( unsigned int i=0; i<numFeatures; i++ )
        {
        sumList[valC][i] += v[i];
        for( unsigned int j=0; j<numFeatures; j++ )
          {
          sumOfSquaresList[valC][i][j] += v[i]*v[j];
          }
        }
      ++countList[valC];
      }

    ++itInMask;
    }

  for( unsigned int c=0; c<numClasses; c++ )
    {
    if( countList[c] > 0 )
      {
      for( unsigned int i=0; i<numFeatures; i++ )
        {
        m_ObjectMeanList[c][i] = sumList[c][i] / countList[c];
        }
      for( unsigned int i=0; i<numFeatures; i++ )
        {
        for( unsigned int j=0; j<numFeatures; j++ )
          {
          m_ObjectCovarianceList[c][i][j] = ( sumOfSquaresList[c][i][j]
            / countList[c] ) - ( m_ObjectMeanList[c][i]
            * m_ObjectMeanList[c][j] );
          }
        }
      }
    }

  timeCollector.Stop( "GenerateStatistics" );

  timeCollector.Report();
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::GenerateLDA()
{
  itk::TimeProbesCollectorBase timeCollector;

  timeCollector.Start( "GenerateLDA" );

  unsigned int numClasses = this->GetNumberOfObjects();
  unsigned int numFeatures = this->GetNumberOfFeatures();

  if( m_PerformLDA )
    {
    ObjectCovarianceType covOfMeans( numFeatures, numFeatures );
    covOfMeans.fill( 0 );
    ObjectCovarianceType meanCov( numFeatures, numFeatures );
    meanCov.fill( 0 );
    for( unsigned int c=0; c<numClasses; c++ )
      {
      for( unsigned int i=0; i<numFeatures; i++ )
        {
        for( unsigned int j=0; j<numFeatures; j++ )
          {
          covOfMeans[i][j] += ( m_ObjectMeanList[c][i] - m_GlobalMean[i] )
            * ( m_ObjectMeanList[c][j] - m_GlobalMean[j] );
          meanCov[i][j] += m_ObjectCovarianceList[c][i][j];
          }
        }
      }

    for( unsigned int i=0; i<numFeatures; i++ )
      {
      for( unsigned int j=0; j<numFeatures; j++ )
        {
        covOfMeans[i][j] /= numClasses;
        meanCov[i][j] /= numClasses;
        }
      }

    ObjectCovarianceType H;
    H = vnl_matrix_inverse<double>(meanCov) * covOfMeans;

    // true = re-order by abs(eval) - zeros will be at end
    ::tube::Eigen<double>( H, m_LDAMatrix, m_LDAValues, true, false );
    }
  else if( m_PerformPCA )
    {
    // true = re-order by abs(eval) - zeros will be at end
    ::tube::Eigen<double>( m_GlobalCovariance, m_LDAMatrix, m_LDAValues,
      true, false );
    }
  else
    {
    ObjectCovarianceType covOfMeans( numFeatures, numFeatures );
    covOfMeans.fill( 0 );
    for( unsigned int c=0; c<numClasses; c++ )
      {
      for( unsigned int i=0; i<numFeatures; i++ )
        {
        for( unsigned int j=0; j<numFeatures; j++ )
          {
          covOfMeans[i][j] += ( m_ObjectMeanList[c][i] - m_GlobalMean[i] )
            * ( m_ObjectMeanList[c][j] - m_GlobalMean[j] );
          }
        }
      }

    for( unsigned int i=0; i<numFeatures; i++ )
      {
      for( unsigned int j=0; j<numFeatures; j++ )
        {
        covOfMeans[i][j] /= numClasses;
        }
      }

    // true = re-order by abs(eval) - zeros will be at end
    ::tube::Eigen<double>( covOfMeans, m_LDAMatrix, m_LDAValues, true,
      false );
    }

  timeCollector.Stop( "GenerateLDA" );

  timeCollector.Report();
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::Update()
{
  this->GenerateStatistics();
  this->GenerateLDA();
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::UpdateLDAImages()
{
}

template <class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_PerformLDA )
    {
    os << indent << "PerformLDA = true" << std::endl;
    }
  else
    {
    os << indent << "PerformLDA = false" << std::endl;
    }
  if( m_PerformPCA )
    {
    os << indent << "PerformPCA = true" << std::endl;
    }
  else
    {
    os << indent << "PerformPCA = false" << std::endl;
    }
  os << indent << "FeatureImageList.size = " << m_FeatureImageList.size()
    << std::endl;
  if( m_Labelmap )
    {
    os << indent << "Labelmap = " << m_Labelmap << std::endl;
    }
  else
    {
    os << indent << "Labelmap = NULL" << std::endl;
    }
  os << indent << "ObjectIdList.size = " << m_ObjectIdList.size()
    << std::endl;
  os << indent << "ObjectMeanList.size = " << m_ObjectMeanList.size()
    << std::endl;
  os << indent << "ObjectCovarianceList.size = "
    << m_ObjectCovarianceList.size() << std::endl;
  os << indent << "GlobalMean = " << m_GlobalMean << std::endl;
  os << indent << "GlobalCovariance = " << m_GlobalCovariance << std::endl;
  os << indent << "LDAMatrix = " << m_LDAMatrix << std::endl;
  os << indent << "LDAValues = " << m_LDAValues << std::endl;
  if( m_LDAImage.IsNotNull() )
    {
    os << indent << "LDAImage = " << m_LDAImage << std::endl;
    }
  else
    {
    os << indent << "LDAImage = NULL" << std::endl;
    }
  if( m_ProgressProcessInfo != NULL )
    {
    os << indent << "ProgressProcessInfo = Set" << std::endl;
    }
  else
    {
    os << indent << "ProgressProcessInfo = NULL" << std::endl;
    }
  os << indent << "ProgressFraction = " << m_ProgressFraction << std::endl;
  os << indent << "ProgressStart = " << m_ProgressStart << std::endl;

}

}

}

#endif //LDAGenerator_txx
