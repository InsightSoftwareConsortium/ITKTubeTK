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

#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkDiscreteGaussianImageFilter.h"

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

  m_LDAUpToDate = false;
  m_LDAImageListUpToDate = false;

  m_FeatureImageList.clear();

  m_Labelmap = NULL;

  m_ObjectIdList.clear();
  m_ObjectIdList.push_back( 1 );
  m_ObjectMeanList.clear();
  m_ObjectCovarianceList.clear();

  m_NumberOfLDA = 0;
  m_LDAMatrix.fill( 0 );
  m_LDAValues.fill( 0 );

  m_LDAImageList.clear();

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
  m_LDAUpToDate = false;
  m_FeatureImageList.clear();
  m_FeatureImageList.push_back( img );
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::AddFeatureImage( typename ImageType::Pointer img )
{
  m_LDAUpToDate = false;
  m_FeatureImageList.push_back( img );
}

template < class ImageT, class LabelmapT >
int
LDAGenerator< ImageT, LabelmapT >
::GetNumberOfFeatureImages( void )
{
  return m_FeatureImageList.size();
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::SetObjectId( ObjectIdType objectId )
{
  m_LDAUpToDate = false;
  m_ObjectIdList.clear();
  m_ObjectIdList.push_back( objectId );
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::AddObjectId( ObjectIdType objectId )
{
  m_LDAUpToDate = false;
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
  if( m_LDAUpToDate )
    {
    return m_NumberOfLDA;
    }
  else
    {
    return 0;
    }
}

template < class ImageT, class LabelmapT >
typename LDAGenerator< ImageT, LabelmapT >::LDAVectorType
LDAGenerator< ImageT, LabelmapT >
::GetLDAVector( unsigned int ldaNum )
{
  if( m_LDAUpToDate && ldaNum < m_NumberOfLDA )
    {
    return m_LDAMatrix.get_column( ldaNum );
    }
  else
    {
    LDAVectorType zero( m_NumberOfLDA );
    zero.fill( 0 );
    return zero;
    }
}

template < class ImageT, class LabelmapT >
double
LDAGenerator< ImageT, LabelmapT >
::GetLDAValue( unsigned int ldaNum )
{
  if( m_LDAUpToDate && ldaNum < m_NumberOfLDA )
    {
    return m_LDAValues[ ldaNum ];
    }
  else
    {
    return 0;
    }
}

template < class ImageT, class LabelmapT >
typename LDAGenerator< ImageT, LabelmapT >::LDAMatrixType *
LDAGenerator< ImageT, LabelmapT >
::GetLDAMatrix( void )
{
  if( m_LDAUpToDate )
    {
    return & m_LDAMatrix;
    }
  else
    {
    return NULL;
    }
}

template < class ImageT, class LabelmapT >
typename LDAGenerator< ImageT, LabelmapT >::LDAValuesType *
LDAGenerator< ImageT, LabelmapT >
::GetLDAValues( void )
{
  if( m_LDAUpToDate )
    {
    return & m_LDAValues;
    }
  else
    {
    return NULL;
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
const typename LDAGenerator< ImageT, LabelmapT >::LDAImageType::Pointer
LDAGenerator< ImageT, LabelmapT >
::GetLDAImage( unsigned int ldaNum )
{
  if( !m_LDAImageListUpToDate || !m_LDAUpToDate )
    {
    if( !m_LDAUpToDate )
      {
      this->Update();
      }
    this->GenerateLDAImages();
    }

  if( ldaNum < m_LDAImageList.size()  )
    {
    return m_LDAImageList[ ldaNum ];
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
void
LDAGenerator< ImageT, LabelmapT >
::GenerateStatistics( void )
{
  itk::TimeProbesCollectorBase timeCollector;

  timeCollector.Start( "GenerateStatistics" );

  typedef itk::ImageRegionConstIterator< MaskImageType >
    ConstMaskImageIteratorType;
  typedef itk::ImageRegionConstIterator< ImageType >
    ConstImageIteratorType;

  ConstMaskImageIteratorType itInMask( m_Labelmap,
    m_Labelmap->GetLargestPossibleRegion() );
  itInMask.GoToBegin();

  unsigned int numClasses = m_ObjectIdList.size();
  unsigned int numFeatures = m_FeatureImageList.size();

  std::vector< ConstImageIteratorType * > itInIm( numFeatures );
  for( unsigned int i=0; i<numFeatures; i++ )
    {
    itInIm[i] = new ConstImageIteratorType( m_FeatureImageList[i],
      m_FeatureImageList[i]->GetLargestPossibleRegion() );
    itInIm[i]->GoToBegin();
    }

  m_ObjectMeanList.resize( numClasses );
  m_ObjectCovarianceList.resize( numClasses );
  std::vector< unsigned int > countList( numClasses );
  ObjectMeanListType sumList( numClasses );
  ObjectCovarianceListType sumOfSquaresList( numClasses );
  for( unsigned int i=0; i<numClasses; i++ )
    {
    m_ObjectMeanList[i].set_size( numFeatures );
    m_ObjectCovarianceList[i].set_size( numFeatures, numFeatures );
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

  FeatureVectorType v( numFeatures );
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
      for( unsigned int i=0; i<numFeatures; i++ )
        {
        v[i] = static_cast< FeatureType >( itInIm[i]->Get() );
        sumList[valC][i] += v[i];
        for( unsigned int j=0; j<numFeatures; j++ )
          {
          sumOfSquaresList[valC][i][j] += v[i]*v[j];
          }
        }
      ++countList[valC];
      }

    ++itInMask;
    for( unsigned int i=0; i<numFeatures; i++ )
      {
      ++( *( itInIm[i] ) );
      }
    }

  for( unsigned int i=0; i<numFeatures; i++ )
    {
    delete itInIm[i];
    }

  for( unsigned int c=0; c<numClasses; c++ )
    {
    globalCount += countList[c];
    for( unsigned int i=0; i<numFeatures; i++ )
      {
      globalSum[i] += sumList[c][i];
      for( unsigned int j=0; j<numFeatures; j++ )
        {
        globalSumOfSquares[i][j] += sumOfSquaresList[c][i][j];
        }
      }
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

  unsigned int numClasses = m_ObjectIdList.size();
  unsigned int numFeatures = m_FeatureImageList.size();

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
  
  m_NumberOfLDA = m_LDAValues.size();
  for( unsigned int i=0; i<m_LDAValues.size(); i++ )
    {
    if( m_LDAValues[i] == 0 )
      {
      m_NumberOfLDA = i;
      break;
      }
    }

  m_LDAUpToDate = true;
  m_LDAImageListUpToDate = false;

  timeCollector.Stop( "GenerateLDA" );

  timeCollector.Report();
}

template < class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::GenerateLDAImages()
{
  if( !m_LDAUpToDate )
    {
    this->Update();
    }

  itk::TimeProbesCollectorBase timeCollector;

  timeCollector.Start( "GenerateLDAImages" );

  typedef itk::ImageRegionConstIterator< ImageType >
    ConstImageIteratorType;
  typedef itk::ImageRegionIterator< LDAImageType >
    ImageIteratorType;

  unsigned int numFeatures = m_FeatureImageList.size();

  std::vector< ConstImageIteratorType * > itInIm( numFeatures );
  for( unsigned int i=0; i<numFeatures; i++ )
    {
    itInIm[i] = new ConstImageIteratorType( m_FeatureImageList[i],
      m_FeatureImageList[i]->GetLargestPossibleRegion() );
    itInIm[i]->GoToBegin();
    }

  typename LDAImageType::RegionType region;
  region = m_FeatureImageList[0]->GetLargestPossibleRegion();
  m_LDAImageList.resize( m_NumberOfLDA );
  for( unsigned int i=0; i<m_NumberOfLDA; i++ )
    {
    m_LDAImageList[i] = LDAImageType::New();
    m_LDAImageList[i]->SetRegions( region );
    m_LDAImageList[i]->CopyInformation( m_FeatureImageList[0] );
    m_LDAImageList[i]->Allocate();
    }

  std::vector< ImageIteratorType * > itLDAIm( m_NumberOfLDA );
  for( unsigned int i=0; i<m_NumberOfLDA; i++ )
    {
    itLDAIm[i] = new ImageIteratorType( m_LDAImageList[i],
      m_LDAImageList[i]->GetLargestPossibleRegion() );
    itLDAIm[i]->GoToBegin();
    }

  FeatureVectorType v( numFeatures );
  FeatureVectorType vLDA( numFeatures );
  while( !itInIm[0]->IsAtEnd() )
    {
    for( unsigned int i=0; i<numFeatures; i++ )
      {
      v[i] = static_cast< FeatureType >( itInIm[i]->Get() );
      ++( *( itInIm[i] ) );
      }

    vLDA = v * m_LDAMatrix; 

    for( unsigned int i=0; i<m_NumberOfLDA; i++ )
      {
      itLDAIm[i]->Set( vLDA[i] );
      ++( *( itLDAIm[i] ) );
      }
    }

  for( unsigned int i=0; i<numFeatures; i++ )
    {
    delete itInIm[i];
    }
  for( unsigned int i=0; i<m_NumberOfLDA; i++ )
    {
    delete itLDAIm[i];
    }

  timeCollector.Stop( "GenerateLDAImages" );

  m_LDAImageListUpToDate = true;

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
  if( !m_LDAUpToDate )
    {
    this->Update();
    }

  this->GenerateLDAImages();
}

template <class ImageT, class LabelmapT >
void
LDAGenerator< ImageT, LabelmapT >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_LDAUpToDate )
    {
    os << indent << "LDAUpToDate = true" << std::endl;
    }
  else
    {
    os << indent << "LDAUpToDate = false" << std::endl;
    }
  if( m_LDAImageListUpToDate )
    {
    os << indent << "LDAImageListUpToDate = true" << std::endl;
    }
  else
    {
    os << indent << "LDAImageListUpToDate = false" << std::endl;
    }

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
  os << indent << "NumberOfLDA = " << m_NumberOfLDA << std::endl;
  os << indent << "LDAMatrix = " << m_LDAMatrix << std::endl;
  os << indent << "LDAValues = " << m_LDAValues << std::endl;
  os << indent << "LDAImageList.size = " << m_LDAImageList.size() 
    << std::endl;
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
