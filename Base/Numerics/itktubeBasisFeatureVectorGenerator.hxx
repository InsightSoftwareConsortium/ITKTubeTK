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

#ifndef __itktubeBasisFeatureVectorGenerator_hxx
#define __itktubeBasisFeatureVectorGenerator_hxx

#include "itktubeBasisFeatureVectorGenerator.h"
#include "tubeMatrixMath.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkTimeProbesCollectorBase.h>

#include <iostream>
#include <limits>

namespace itk
{

namespace tube
{

template< class TImage, class TLabelMap >
BasisFeatureVectorGenerator< TImage, TLabelMap >
::BasisFeatureVectorGenerator( void )
{
  m_LabelMap = NULL;

  m_ObjectIdList.clear();
  m_ObjectMeanList.clear();
  m_ObjectCovarianceList.clear();

  m_GlobalMean = 0;
  m_GlobalCovariance = 0;

  m_NumberOfBasis = 0;

  m_NumberOfBasisToUseAsFeatures = 0;
  m_NumberOfLDABasisToUseAsFeatures = 0;

  m_InputFeatureVectorGenerator = NULL;

  m_BasisValues.set_size( 0 );
  m_BasisMatrix.set_size( 0, 0 );
}

template< class TImage, class TLabelMap >
BasisFeatureVectorGenerator< TImage, TLabelMap >
::~BasisFeatureVectorGenerator( void )
{
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::SetInputFeatureVectorGenerator( FeatureVectorGeneratorType * fGen )
{
  m_InputFeatureVectorGenerator = fGen;
}

template< class TImage, class TLabelMap >
typename FeatureVectorGenerator< TImage >::Pointer
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetInputFeatureVectorGenerator( void )
{
  return m_InputFeatureVectorGenerator;
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::SetObjectId( ObjectIdType objectId )
{
  m_ObjectIdList.clear();
  m_ObjectIdList.push_back( objectId );
  m_NumberOfLDABasisToUseAsFeatures = 1;
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::AddObjectId( ObjectIdType objectId )
{
  m_ObjectIdList.push_back( objectId );
  ++m_NumberOfLDABasisToUseAsFeatures;
}

template< class TImage, class TLabelMap >
typename BasisFeatureVectorGenerator< TImage, TLabelMap >::ObjectIdType
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetObjectId( unsigned int num ) const
{
  if( num < m_ObjectIdList.size() )
    {
    return m_ObjectIdList[num];
    }
  else
    {
    throw;
    }
}

template< class TImage, class TLabelMap >
unsigned int
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetNumberOfObjectIds( void ) const
{
  return m_ObjectIdList.size();
}

template< class TImage, class TLabelMap >
typename BasisFeatureVectorGenerator< TImage, TLabelMap >::ValueType
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetObjectMean( ObjectIdType num ) const
{
  if( num < m_ObjectIdList.size() )
    {
    return m_ObjectMeanList[num];
    }
  else
    {
    throw;
    }
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::SetObjectMean( ObjectIdType num, ValueType val )
{
  if( num < m_ObjectIdList.size() )
    {
    m_ObjectMeanList[num] = val;
    }
}

template< class TImage, class TLabelMap >
typename BasisFeatureVectorGenerator< TImage, TLabelMap >::MatrixType
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetObjectCovariance( ObjectIdType num ) const
{
  if( num < m_ObjectIdList.size() )
    {
    return m_ObjectCovarianceList[num];
    }
  else
    {
    throw;
    }
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::SetObjectCovariance( ObjectIdType num, MatrixType val )
{
  if( num < m_ObjectIdList.size() )
    {
    m_ObjectCovarianceList[num] = val;
    }
  else
    {
    throw;
    }
}

template< class TImage, class TLabelMap >
typename BasisFeatureVectorGenerator< TImage, TLabelMap >::ValueType
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetGlobalMean( void ) const
{
  return m_GlobalMean;
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::SetGlobalMean( ValueType val )
{
  m_GlobalMean = val;
}

template< class TImage, class TLabelMap >
typename BasisFeatureVectorGenerator< TImage, TLabelMap >::MatrixType
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetGlobalCovariance( void ) const
{
  return m_GlobalCovariance;
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::SetGlobalCovariance( MatrixType val )
{
  m_GlobalCovariance = val;
}


template< class TImage, class TLabelMap >
unsigned int
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetNumberOfBasis( void ) const
{
  return m_BasisValues.size();
}

template< class TImage, class TLabelMap >
typename BasisFeatureVectorGenerator< TImage, TLabelMap >::VectorType
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetBasisVector( unsigned int basisNum ) const
{
  if( basisNum < m_BasisValues.size() )
    {
    return m_BasisMatrix.get_column( basisNum );
    }
  else
    {
    std::cerr << "Basis does not exist." << std::endl;
    return m_BasisMatrix.get_column( 0 );
    }
}

template< class TImage, class TLabelMap >
double
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetBasisValue( unsigned int basisNum ) const
{
  if( basisNum < m_BasisValues.size() )
    {
    return m_BasisValues[basisNum];
    }
  else
    {
    std::cerr << "Basis does not exist." << std::endl;
    return m_BasisValues[basisNum];
    }
}

template< class TImage, class TLabelMap >
typename BasisFeatureVectorGenerator< TImage, TLabelMap >::MatrixType
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetBasisMatrix( void ) const
{
  return m_BasisMatrix;
}

template< class TImage, class TLabelMap >
typename BasisFeatureVectorGenerator< TImage, TLabelMap >::VectorType
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetBasisValues( void ) const
{
  return m_BasisValues;
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::SetBasisValue( unsigned int basisNum, double value )
{
  if( basisNum < m_BasisValues.size() )
    {
    m_BasisValues[basisNum] = value;
    }
  else
    {
    throw;
    }
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::SetBasisMatrix( const MatrixType & mat )
{
  m_BasisMatrix = mat;
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
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

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::SetBasisValues( const VectorType & values )
{
  m_BasisValues = values;
}

template< class TImage, class TLabelMap >
typename BasisFeatureVectorGenerator< TImage, TLabelMap >::FeatureImageType::Pointer
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetFeatureImage( unsigned int featureNum ) const
{
  if( featureNum < m_NumberOfBasisToUseAsFeatures  )
    {
    itk::TimeProbesCollectorBase timeCollector;

    timeCollector.Start( "GenerateBasisImage" );

    typedef itk::ImageRegionIteratorWithIndex< FeatureImageType >
      ImageIteratorType;

    typename FeatureImageType::RegionType region;
    region = this->m_InputImageList[0]->GetLargestPossibleRegion();

    typename FeatureImageType::Pointer featureImage =
      FeatureImageType::New();
    featureImage->SetRegions( region );
    featureImage->CopyInformation( this->m_InputImageList[0] );
    featureImage->Allocate();

    ImageIteratorType itBasisIm( featureImage,
      featureImage->GetLargestPossibleRegion() );

    if( m_LabelMap.IsNotNull() )
      {
      const unsigned int numClasses = this->GetNumberOfObjectIds();
      typedef itk::ImageRegionConstIteratorWithIndex< LabelMapType >
        ConstLabelMapIteratorType;
      ConstLabelMapIteratorType itInMask( m_LabelMap,
        m_LabelMap->GetLargestPossibleRegion() );
      bool found = false;
      ObjectIdType previousMaskValue =
        static_cast<ObjectIdType>( itInMask.Get() ) + 1;
      while( !itBasisIm.IsAtEnd() )
        {
        ObjectIdType maskVal = static_cast<ObjectIdType>( itInMask.Get() );
        if( maskVal != previousMaskValue )
          {
          found = false;
          previousMaskValue = maskVal;
          for( unsigned int c = 0; c < numClasses; c++ )
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
          itBasisIm.Set( this->GetFeatureVectorValue( indx, featureNum ) );
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
        itBasisIm.Set( this->GetFeatureVectorValue( indx, featureNum ) );
        ++itBasisIm;
        }
      }

    timeCollector.Stop( "GenerateBasisImage" );
    timeCollector.Report();

    return featureImage;
    }
  else
    {
    throw;
    }
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GenerateBasis( void )
{
  itk::TimeProbesCollectorBase timeCollector;

  timeCollector.Start( "GenerateStatistics" );

  typedef itk::ImageRegionConstIteratorWithIndex< LabelMapType >
    ConstLabelMapIteratorType;
  ConstLabelMapIteratorType itInMask( m_LabelMap,
    m_LabelMap->GetLargestPossibleRegion() );

  const unsigned int numClasses = this->GetNumberOfObjectIds();
  const unsigned int numInputFeatures =
    m_InputFeatureVectorGenerator->GetNumberOfFeatures();
  m_ObjectMeanList.resize( numClasses );
  m_ObjectCovarianceList.resize( numClasses );
  std::vector< unsigned int > countList( numClasses );
  for( unsigned int i = 0; i < numClasses; i++ )
    {
    m_ObjectMeanList[i].set_size( numInputFeatures );
    m_ObjectMeanList[i].fill( 0 );
    m_ObjectCovarianceList[i].set_size( numInputFeatures,
      numInputFeatures );
    m_ObjectCovarianceList[i].fill( 0 );
    countList[i] = 0;
    }

  m_GlobalMean.set_size( numInputFeatures );
  m_GlobalMean.fill( 0 );
  m_GlobalCovariance.set_size( numInputFeatures, numInputFeatures );
  m_GlobalCovariance.fill( 0 );

  unsigned int globalCount = 0;

  double delta = 0;
  double deltaJ = 0;
  unsigned int valC = 0;
  bool found = false;
  itInMask.GoToBegin();
  ObjectIdType previousMaskValue =
    static_cast<ObjectIdType>( itInMask.Get() ) + 1;
  while( !itInMask.IsAtEnd() )
    {
    ObjectIdType maskVal = static_cast<ObjectIdType>( itInMask.Get() );
    if( maskVal != previousMaskValue )
      {
      previousMaskValue = maskVal;
      found = false;
      for( unsigned int c = 0; c < numClasses; c++ )
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
      FeatureVectorType v = m_InputFeatureVectorGenerator->
        GetFeatureVector( indx );

      ++globalCount;
      ++countList[valC];
      for( unsigned int i = 0; i < numInputFeatures; i++ )
        {
        delta = v[i] - m_GlobalMean[i];
        m_GlobalMean[i] += delta / globalCount;

        delta = v[i] - m_ObjectMeanList[valC][i];
        m_ObjectMeanList[valC][i] += delta / countList[valC];
        }
      for( unsigned int i = 0; i < numInputFeatures; i++ )
        {
        for( unsigned int j = 0; j < numInputFeatures; j++ )
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

  for( unsigned int i = 0; i < numInputFeatures; i++ )
    {
    for( unsigned int j = 0; j < numInputFeatures; j++ )
      {
      m_GlobalCovariance[i][j] /= globalCount - 1;
      for( unsigned int c = 0; c < numClasses; c++ )
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

  // PCA
  ::tube::ComputeEigen<double>( m_GlobalCovariance,
    m_BasisMatrix, m_BasisValues, true, false );

  if( m_NumberOfLDABasisToUseAsFeatures > 0 )
    {
    MatrixType covarianceOfMeans( numInputFeatures, numInputFeatures );
    covarianceOfMeans.fill( 0 );
    MatrixType meanCovariance( numInputFeatures, numInputFeatures );
    meanCovariance.fill( 0 );
    for( unsigned int c = 0; c < numClasses; c++ )
      {
      for( unsigned int i = 0; i < numInputFeatures; i++ )
        {
        for( unsigned int j = 0; j < numInputFeatures; j++ )
          {
          covarianceOfMeans[i][j] +=
            ( m_ObjectMeanList[c][i] - m_GlobalMean[i] )
            * ( m_ObjectMeanList[c][j] - m_GlobalMean[j] );
          meanCovariance[i][j] += ( m_ObjectCovarianceList[c][i][j]
            * m_ObjectCovarianceList[c][i][j] );
          }
        }
      }

    for( unsigned int i = 0; i < numInputFeatures; i++ )
      {
      for( unsigned int j = 0; j < numInputFeatures; j++ )
        {
        covarianceOfMeans[i][j] /= numClasses;
        meanCovariance[i][j] = vcl_sqrt( meanCovariance[i][j] / numClasses );
        }
      }

    MatrixType H;
    H = covarianceOfMeans * vnl_matrix_inverse<double>( meanCovariance );

    MatrixType ldaBasisMatrix;
    VectorType ldaBasisValues;
    ::tube::ComputeEigen<double>( H, ldaBasisMatrix, ldaBasisValues,
      true, false );
    for( int f = m_BasisValues.size()-1;
      f >= (int)m_NumberOfLDABasisToUseAsFeatures; --f )
      {
      m_BasisValues[f] = m_BasisValues[ f - m_NumberOfLDABasisToUseAsFeatures ];
      m_BasisMatrix.set_column( f, m_BasisMatrix.get_column( f -
        m_NumberOfLDABasisToUseAsFeatures ) );
      }
    for( unsigned int f = 0; f < m_NumberOfLDABasisToUseAsFeatures; ++f )
      {
      m_BasisValues[f] = ldaBasisValues[ f ];
      m_BasisMatrix.set_column( f, ldaBasisMatrix.get_column( f  ) );
      }
    }

  m_NumberOfBasisToUseAsFeatures = m_BasisValues.size();

  timeCollector.Stop( "GenerateBasis" );

  timeCollector.Report();
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::SetNumberOfBasisToUseAsFeatures( unsigned int numBasisUsed )
{
  m_NumberOfBasisToUseAsFeatures = numBasisUsed;
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::SetNumberOfLDABasisToUseAsFeatures( unsigned int numBasisUsed )
{
  m_NumberOfLDABasisToUseAsFeatures = numBasisUsed;
}

template< class TImage, class TLabelMap >
unsigned int
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetNumberOfFeatures( void ) const
{
  return m_NumberOfBasisToUseAsFeatures;
}

template< class TImage, class TLabelMap >
typename BasisFeatureVectorGenerator< TImage, TLabelMap >::FeatureVectorType
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetFeatureVector( const IndexType & indx ) const
{
  const unsigned int numFeatures = this->GetNumberOfFeatures();

  FeatureVectorType featureVector;
  featureVector.set_size( numFeatures );

  VectorType vBasis;
  FeatureVectorType vInput;

  for( unsigned int i = 0; i < numFeatures; i++ )
    {
    vBasis = this->GetBasisVector( i );
    vInput = m_InputFeatureVectorGenerator->GetFeatureVector( indx );
    featureVector[i] = 0;
    for( unsigned int j = 0; j < vBasis.size(); j++ )
      {
      featureVector[i] += vBasis[j] * vInput[j];
      }
    }

  return featureVector;
}

template< class TImage, class TLabelMap >
typename BasisFeatureVectorGenerator< TImage, TLabelMap >::FeatureValueType
BasisFeatureVectorGenerator< TImage, TLabelMap >
::GetFeatureVectorValue( const IndexType & indx,
  unsigned int featureNum ) const
{
  VectorType vBasis;
  FeatureVectorType vInput;

  vBasis = this->GetBasisVector( featureNum );
  vInput = m_InputFeatureVectorGenerator->GetFeatureVector( indx );

  FeatureValueType featureVector = 0;
  for( unsigned int j = 0; j < vBasis.size(); j++ )
    {
    featureVector += vBasis[j] * vInput[j];
    }
  return featureVector;
}

template< class TImage, class TLabelMap >
void
BasisFeatureVectorGenerator< TImage, TLabelMap >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_LabelMap )
    {
    os << indent << "LabelMap = " << m_LabelMap << std::endl;
    }
  else
    {
    os << indent << "LabelMap = NULL" << std::endl;
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

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeBasisFeatureVectorGenerator_hxx)
