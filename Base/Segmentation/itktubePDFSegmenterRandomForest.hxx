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

#ifndef __itktubePDFSegmenterRandomForest_hxx
#define __itktubePDFSegmenterRandomForest_hxx

#include "itktubePDFSegmenterRandomForest.h"

namespace itk
{

namespace tube
{

template< class TImage, class TLabelMap >
PDFSegmenterRandomForest< TImage, TLabelMap >
::PDFSegmenterRandomForest( void )
{
  m_TrainingDataStride = 1;

  m_NumberOfDecisionTrees = 100;
}

template< class TImage, class TLabelMap >
PDFSegmenterRandomForest< TImage, TLabelMap >
::~PDFSegmenterRandomForest( void )
{
}

template< class TImage, class TLabelMap >
typename PDFSegmenterRandomForest< TImage, TLabelMap>::DecisionForestType &
PDFSegmenterRandomForest< TImage, TLabelMap >
::GetModel( void )
{
  return m_Model;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterRandomForest< TImage, TLabelMap >
::SetModel( typename PDFSegmenterRandomForest< TImage, TLabelMap>::
  DecisionForestType & model )
{
  m_Model = model;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterRandomForest< TImage, TLabelMap >
::GeneratePDFs( void )
{
  if( !this->m_SampleUpToDate )
    {
    this->GenerateSample();
    }
  this->m_PDFsUpToDate = true;

  unsigned int numClasses = this->m_ObjectIdList.size();
  unsigned int numFeatures = this->m_FeatureVectorGenerator->
    GetNumberOfFeatures();

  unsigned int sampleSize = 0;
  for( unsigned int c=0; c<numClasses; ++c )
    {
    sampleSize += this->m_InClassList[c].size();
    }

  const size_t shape[] = { sampleSize, numFeatures };
  andres::Marray< FeatureValueType > features( shape, shape + 2 );
  andres::Marray< unsigned int > labels( shape, shape + 1 );

  unsigned int sampleNum = 0;
  for( unsigned int c=0; c<numClasses; ++c )
    {
    typename ListSampleType::const_iterator
      inClassListIt( this->m_InClassList[c].begin() );
    typename ListSampleType::const_iterator
      inClassListItEnd( this->m_InClassList[c].end() );
    while( sampleNum < sampleSize && inClassListIt != inClassListItEnd )
      {
      for( unsigned int f=0; f<numFeatures; ++f )
        {
        features( sampleNum, f ) = (*inClassListIt)[ f ];
        }
      labels( sampleNum ) = c;
      ++sampleNum;
      for( unsigned int s = 0; inClassListIt != inClassListItEnd
        && s < m_TrainingDataStride; ++s )
        {
        ++inClassListIt;
        }
      }
    }

  std::cout << "Learn" << std::endl;
  m_Model.learn( features, labels, m_NumberOfDecisionTrees );
  std::cout << "...done" << std::endl;
}

template< class TImage, class TLabelMap >
typename PDFSegmenterRandomForest< TImage, TLabelMap >
::ProbabilityVectorType
PDFSegmenterRandomForest< TImage, TLabelMap >
::GetProbabilityVector( const FeatureVectorType & fv) const
{
  unsigned int numClasses = this->m_ObjectIdList.size();
  unsigned int numFeatures = this->GetNumberOfFeatures();

  const size_t shape[] = { 1, numFeatures };
  andres::Marray< FeatureValueType > features( shape, shape + 2 );

  for( unsigned int f=0; f<numFeatures; ++f )
    {
    features( 0, f ) = fv[ f ];
    }

  andres::Marray< ProbabilityPixelType > probabilities( shape, shape + 2 );

  m_Model.predict( features, probabilities );

  ProbabilityVectorType prob( numClasses );
  for( unsigned int c=0; c<numClasses; ++c )
    {
    prob[c] = probabilities( c );
    }

  return prob;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterRandomForest< TImage, TLabelMap >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Training data stride = "
    << m_TrainingDataStride << std::endl;

  os << indent << "Number of decision trees = "
    << m_NumberOfDecisionTrees << std::endl;

  os << indent << "Model = model" << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubePDFSegmenterRandomForest_hxx)
