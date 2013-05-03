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
#ifndef __itkTubeBasisFeatureVectorGenerator_txx
#define __itkTubeBasisFeatureVectorGenerator_txx

#include <limits>
#include <iostream>

#include "itkTubeBasisFeatureVectorGenerator.h"

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
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::BasisFeatureVectorGenerator()
{
  m_PerformLDA = true;
  m_PerformPCA = false;

  m_Labelmap = NULL;

  m_ObjectIdList.clear();
  m_ObjectMeanList.clear();
  m_ObjectCovarianceList.clear();

  m_GlobalMean = 0;
  m_GlobalCovariance = 0;

  m_NumberOfBasisToUseAsFeatures = 0;
  m_InputFeatureVectorGenerator = NULL;

  m_BasisValues.set_size( 0 );
  m_BasisMatrix.set_size( 0, 0 );
}

template< class ImageT, class LabelmapT >
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::~BasisFeatureVectorGenerator()
{
}

template < class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::SetInputFeatureVectorGenerator( FeatureVectorGeneratorType * fGen )
{
  m_InputFeatureVectorGenerator = fGen;
}

template < class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::SetObjectId( ObjectIdType objectId )
{
  m_ObjectIdList.clear();
  m_ObjectIdList.push_back( objectId );
}

template < class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::AddObjectId( ObjectIdType objectId )
{
  m_ObjectIdList.push_back( objectId );
}

template < class ImageT, class LabelmapT >
typename BasisFeatureVectorGenerator< ImageT, LabelmapT >::ObjectIdType
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetObjectId( unsigned int num ) const
{
  if( num < m_ObjectIdList.size() )
    {
    return m_ObjectIdList[ num ];
    }
  else
    {
    throw;
    }
}

template < class ImageT, class LabelmapT >
unsigned int
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetNumberOfObjectIds( void ) const
{
  return m_ObjectIdList.size();
}

template < class ImageT, class LabelmapT >
typename BasisFeatureVectorGenerator< ImageT, LabelmapT >::ValueType
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetObjectMean( ObjectIdType num ) const
{
  if( num < m_ObjectIdList.size() )
    {
    return m_ObjectMeanList[ num ];
    }
  else
    {
    throw;
    }
}

template < class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::SetObjectMean( ObjectIdType num, ValueType val )
{
  if( num < m_ObjectIdList.size() )
    {
    m_ObjectMeanList[ num ] = val;
    }
}

template < class ImageT, class LabelmapT >
typename BasisFeatureVectorGenerator< ImageT, LabelmapT >::MatrixType
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetObjectCovariance( ObjectIdType num ) const
{
  if( num < m_ObjectIdList.size() )
    {
    return m_ObjectCovarianceList[ num ];
    }
  else
    {
    throw;
    }
}

template < class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::SetObjectCovariance( ObjectIdType num, MatrixType val )
{
  if( num < m_ObjectIdList.size() )
    {
    m_ObjectCovarianceList[ num ] = val;
    }
  else
    {
    throw;
    }
}

template < class ImageT, class LabelmapT >
typename BasisFeatureVectorGenerator< ImageT, LabelmapT >::ValueType
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetGlobalMean( void ) const
{
  return m_GlobalMean;
}

template < class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::SetGlobalMean( ValueType val )
{
  m_GlobalMean = val;
}

template < class ImageT, class LabelmapT >
typename BasisFeatureVectorGenerator< ImageT, LabelmapT >::MatrixType
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetGlobalCovariance( void ) const
{
  return m_GlobalCovariance;
}

template < class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::SetGlobalCovariance( MatrixType val )
{
  m_GlobalCovariance = val;
}


template < class ImageT, class LabelmapT >
unsigned int
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetNumberOfBasis( void ) const
{
  return m_BasisValues.size();
}

template < class ImageT, class LabelmapT >
typename BasisFeatureVectorGenerator< ImageT, LabelmapT >::VectorType
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetBasisVector( unsigned int basisNum ) const
{
  if( basisNum < m_BasisValues.size() )
    {
    return m_BasisMatrix.get_column( basisNum );
    }
  else
    {
    throw;
    }
}

template < class ImageT, class LabelmapT >
double
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetBasisValue( unsigned int basisNum ) const
{
  if( basisNum < m_BasisValues.size() )
    {
    return m_BasisValues[ basisNum ];
    }
  else
    {
    throw;
    }
}

template < class ImageT, class LabelmapT >
typename BasisFeatureVectorGenerator< ImageT, LabelmapT >::MatrixType
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetBasisMatrix( void ) const
{
  return m_BasisMatrix;
}

template < class ImageT, class LabelmapT >
typename BasisFeatureVectorGenerator< ImageT, LabelmapT >::VectorType
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetBasisValues( void ) const
{
  return m_BasisValues;
}

template < class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::SetBasisValue( unsigned int basisNum, double value )
{
  if( basisNum < m_BasisValues.size() )
    {
    m_BasisValues[ basisNum ] = value;
    }
  else
    {
    throw;
    }
}

template < class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::SetBasisMatrix( const MatrixType & mat )
{
  m_BasisMatrix = mat;
}

template < class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::SetBasisVector( unsigned int basisNum, const VectorType & vec )
{
  if( basisNum < m_BasisValues.size() )
    {
    m_BasisMatrix.set_column( basisNum, vec );
    }
  else
    {
    throw;
    }
}

template < class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::SetBasisValues( const VectorType & values )
{
  m_BasisValues = values;
}

template < class ImageT, class LabelmapT >
typename BasisFeatureVectorGenerator< ImageT, LabelmapT >::FeatureImageType::Pointer
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetFeatureImage( unsigned int fNum ) const
{
  if( fNum < m_NumberOfBasisToUseAsFeatures  )
    {
    itk::TimeProbesCollectorBase timeCollector;

    timeCollector.Start( "GenerateBasisImage" );

    typedef itk::ImageRegionIteratorWithIndex< FeatureImageType >
      ImageIteratorType;

    unsigned int numFeatures = this->GetNumberOfFeatures();

    typename FeatureImageType::RegionType region;
    region = this->m_InputImageList[0]->GetLargestPossibleRegion();

    typename FeatureImageType::Pointer fi = FeatureImageType::New();
    fi->SetRegions( region );
    fi->CopyInformation( this->m_InputImageList[0] );
    fi->Allocate();

    ImageIteratorType itBasisIm( fi, fi->GetLargestPossibleRegion() );

    FeatureVectorType fv( numFeatures );
    VectorType v( numFeatures );
    VectorType vBasis( numFeatures );

    if( m_Labelmap.IsNotNull() )
      {
      unsigned int numClasses = this->GetNumberOfObjectIds();
      typedef itk::ImageRegionConstIteratorWithIndex< LabelmapType >
        ConstLabelmapIteratorType;
      ConstLabelmapIteratorType itInMask( m_Labelmap,
        m_Labelmap->GetLargestPossibleRegion() );
      bool found = false;
      ObjectIdType prevMaskVal = static_cast<ObjectIdType>(itInMask.Get())
        + 1;
      while( !itBasisIm.IsAtEnd() )
        {
        ObjectIdType maskVal = static_cast<ObjectIdType>( itInMask.Get() );
        if( maskVal != prevMaskVal )
          {
          found = false;
          prevMaskVal = maskVal;
          for( unsigned int c=0; c<numClasses; c++ )
            {
            if( maskVal == m_ObjectIdList[c] )
              {
              found = true;
              break;
              }
            }
          }
        if( !found )
          {
          itBasisIm.Set( 0 );
          }
        else
          {
          IndexType indx = itBasisIm.GetIndex();
          fv = this->GetFeatureVector( indx );
          for( unsigned int f=0; f<numFeatures; f++ )
            {
            v[f] = fv[f];
            }

          vBasis = v * m_BasisMatrix;

          itBasisIm.Set( vBasis[ fNum ] );
          }
        ++itBasisIm;
        ++itInMask;
        }
      }
    else
      {
      while( !itBasisIm.IsAtEnd() )
        {
        IndexType indx = itBasisIm.GetIndex();
        fv = this->GetFeatureVector( indx );
        for( unsigned int f=0; f<numFeatures; f++ )
          {
          v[f] = fv[f];
          }

        vBasis = v * m_BasisMatrix;

        itBasisIm.Set( vBasis[ fNum ] );
        ++itBasisIm;
        }
      }

    timeCollector.Stop( "GenerateBasisImage" );
    timeCollector.Report();

    return fi;
    }
  else
    {
    throw;
    }
}

template < class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GenerateBasis( void )
{
  itk::TimeProbesCollectorBase timeCollector;

  timeCollector.Start( "GenerateStatistics" );

  typedef itk::ImageRegionConstIteratorWithIndex< LabelmapType >
    ConstLabelmapIteratorType;
  ConstLabelmapIteratorType itInMask( m_Labelmap,
    m_Labelmap->GetLargestPossibleRegion() );

  unsigned int numClasses = this->GetNumberOfObjectIds();
  unsigned int numFeatures = this->GetNumberOfFeatures();

  m_ObjectMeanList.resize( numClasses );
  m_ObjectCovarianceList.resize( numClasses );
  std::vector< unsigned int > countList( numClasses );
  for( unsigned int i=0; i<numClasses; i++ )
    {
    m_ObjectMeanList[i].set_size( numFeatures );
    m_ObjectMeanList[i].fill( 0 );
    m_ObjectCovarianceList[i].set_size( numFeatures, numFeatures );
    m_ObjectCovarianceList[i].fill( 0 );
    countList[i] = 0;
    }

  m_GlobalMean.set_size( numFeatures );
  m_GlobalMean.fill( 0 );
  m_GlobalCovariance.set_size( numFeatures, numFeatures );
  m_GlobalCovariance.fill( 0 );

  unsigned int globalCount = 0;

  double delta = 0;
  double deltaJ = 0;
  unsigned int valC = 0;
  bool found = false;
  itInMask.GoToBegin();
  ObjectIdType prevMaskVal = static_cast<ObjectIdType>(itInMask.Get()) + 1;
  while( !itInMask.IsAtEnd() )
    {
    ObjectIdType maskVal = static_cast<ObjectIdType>( itInMask.Get() );
    if( maskVal != prevMaskVal )
      {
      prevMaskVal = maskVal;
      found = false;
      for( unsigned int c=0; c<numClasses; c++ )
        {
        if( maskVal == m_ObjectIdList[c] )
          {
          valC = c;
          found = true;
          break;
          }
        }
      }
    if( found )
      {
      IndexType indx = itInMask.GetIndex();
      FeatureVectorType v = this->GetFeatureVector( indx );

      ++globalCount;
      ++countList[valC];
      for( unsigned int i=0; i<numFeatures; i++ )
        {
        delta = v[i] - m_GlobalMean[i];
        m_GlobalMean[i] += delta / globalCount;

        delta = v[i] - m_ObjectMeanList[valC][i];
        m_ObjectMeanList[valC][i] += delta / countList[valC];
        }
      for( unsigned int i=0; i<numFeatures; i++ )
        {
        for( unsigned int j=0; j<numFeatures; j++ )
          {
          delta = v[i] - m_GlobalMean[i];
          deltaJ = v[j] - m_GlobalMean[j];
          m_GlobalCovariance[i][j] += delta * deltaJ;

          delta = v[i] - m_ObjectMeanList[valC][i];
          deltaJ = v[j] - m_ObjectMeanList[valC][j];
          m_ObjectCovarianceList[valC][i][j] += delta * deltaJ;
          }
        }
      }
    ++itInMask;
    }

  for( unsigned int i=0; i<numFeatures; i++ )
    {
    for( unsigned int j=0; j<numFeatures; j++ )
      {
      m_GlobalCovariance[i][j] /= globalCount - 1;
      for( unsigned int c=0; c<numClasses; c++ )
        {
        if( countList[c] > 1 )
          {
          m_ObjectCovarianceList[c][i][j] /= countList[c] - 1;
          }
        else
          {
          m_ObjectCovarianceList[c][i][j] = 1;
          }
        }
      }
    }

  timeCollector.Stop( "GenerateStatistics" );

  timeCollector.Start( "GenerateBasis" );

  if( m_PerformLDA )
    {
    std::cout << "LDA Global Means = " << m_GlobalMean << std::endl;
    std::cout << "LDA Global Cov = " << m_GlobalCovariance << std::endl;
    MatrixType covOfMeans( numFeatures, numFeatures );
    covOfMeans.fill( 0 );
    MatrixType meanCov( numFeatures, numFeatures );
    meanCov.fill( 0 );
    for( unsigned int c=0; c<numClasses; c++ )
      {
      std::cout << "Object Means[" << m_ObjectIdList[c] << "] = "
        << m_ObjectMeanList[c] << std::endl;
      std::cout << "LDA Cov[" << m_ObjectIdList[c] << "] = "
        << m_ObjectCovarianceList[c] << std::endl;
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

    MatrixType H;
    H = covOfMeans * vnl_matrix_inverse<double>(meanCov);

    // true = re-order by abs(eval) - zeros will be at end
    ::tube::ComputeEigen<double>( H, m_BasisMatrix, m_BasisValues,
      true, false );
    }
  else if( m_PerformPCA )
    {
    // true = re-order by abs(eval) - zeros will be at end
    ::tube::ComputeEigen<double>( m_GlobalCovariance,
      m_BasisMatrix, m_BasisValues, true, false );
    }
  else
    {
    MatrixType covOfMeans( numFeatures, numFeatures );
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
    ::tube::ComputeEigen<double>( covOfMeans, m_BasisMatrix, m_BasisValues,
      true, false );
    }

  timeCollector.Stop( "GenerateBasis" );

  timeCollector.Report();
}

template < class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::SetNumberOfBasisToUseAsFeatures( unsigned int numBasisUsed )
{
  m_NumberOfBasisToUseAsFeatures = numBasisUsed;
}

template < class ImageT, class LabelmapT >
unsigned int
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetNumberOfFeatures( void ) const
{
  return m_NumberOfBasisToUseAsFeatures;
}

template < class ImageT, class LabelmapT >
typename BasisFeatureVectorGenerator< ImageT, LabelmapT >::FeatureVectorType
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetFeatureVector( const IndexType & indx ) const
{
  unsigned int numFeatures = this->GetNumberOfFeatures();

  FeatureVectorType fv;
  fv.set_size( numFeatures );

  VectorType vBasis;
  FeatureVectorType vInput;

  for( unsigned int i=0; i<numFeatures; i++ )
    {
    vBasis = this->GetBasisVector( i );
    vInput = m_InputFeatureVectorGenerator->GetFeatureVector( indx );
    fv[i] = 0;
    for( unsigned int j=0; j<vBasis.size(); j++ )
      {
      fv[i] += vBasis[j] * vInput[j];
      }
    }

  return fv;
}

template < class ImageT, class LabelmapT >
typename BasisFeatureVectorGenerator< ImageT, LabelmapT >::FeatureValueType
BasisFeatureVectorGenerator< ImageT, LabelmapT >
::GetFeatureVectorValue( const IndexType & indx, unsigned int fNum ) const
{
  VectorType vBasis;
  FeatureVectorType vInput;

  vBasis = this->GetBasisVector( fNum );
  vInput = m_InputFeatureVectorGenerator->GetFeatureVector( indx );

  FeatureValueType fv = 0;
  for( unsigned int j=0; j<vBasis.size(); j++ )
    {
    fv += vBasis[j] * vInput[j];
    }
  return fv;
}

template <class ImageT, class LabelmapT >
void
BasisFeatureVectorGenerator< ImageT, LabelmapT >
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
  os << indent << "BasisMatrix = " << m_BasisMatrix << std::endl;
  os << indent << "BasisValues = " << m_BasisValues << std::endl;
}

}

}

#endif //BasisFeatureVectorGenerator_txx
