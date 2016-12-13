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

#ifndef __itktubePDFSegmenterSVM_hxx
#define __itktubePDFSegmenterSVM_hxx

#include "itktubePDFSegmenterSVM.h"

#include "svm.h"

namespace itk
{

namespace tube
{

template< class TImage, class TLabelMap >
PDFSegmenterSVM< TImage, TLabelMap >
::PDFSegmenterSVM( void )
{
  m_TrainingDataStride = 1.0;

  m_Model = NULL;

  m_SVMClassWeight.clear();

  m_Parameter.svm_type = C_SVC;
  m_Parameter.kernel_type = LINEAR; // default: RBF;
  m_Parameter.degree = 3;
  m_Parameter.gamma = 0; // Set to 1.0/num_features in GeneratePDFs
  m_Parameter.coef0 = 0;
  m_Parameter.nu = 0.5;
  m_Parameter.cache_size = 100;
  m_Parameter.C = 1;
  m_Parameter.eps = 0.1; // default: 1e-3;
  m_Parameter.p = 0.1;
  m_Parameter.shrinking = 0;
  m_Parameter.probability = 1;
  m_Parameter.nr_weight = 0;
  m_Parameter.weight_label = NULL;
  m_Parameter.weight = NULL;

  m_Problem.l = 0;
  m_Problem.y = NULL;
  m_Problem.x = NULL;

  m_Space = NULL;
}

template< class TImage, class TLabelMap >
PDFSegmenterSVM< TImage, TLabelMap >
::~PDFSegmenterSVM( void )
{
  if( m_Problem.y != NULL )
    {
    free( m_Problem.y );
    m_Problem.y = NULL;
    }
  if( m_Problem.x != NULL )
    {
    free( m_Problem.x );
    m_Problem.x = NULL;
    }
  if( m_Space != NULL )
    {
    free( m_Space );
    m_Space = NULL;
    }
  if( m_Model != NULL )
    {
    svm_free_and_destroy_model( & m_Model );
    }
  svm_destroy_param( & m_Parameter );
}

template< class TImage, class TLabelMap >
void
PDFSegmenterSVM< TImage, TLabelMap >
::SetSVMClassWeight( unsigned int c, double w )
{
  unsigned int numClasses = this->m_ObjectIdList.size();

  if( m_SVMClassWeight.size() != numClasses )
    {
    m_SVMClassWeight.resize( numClasses );
    for( unsigned int i=0; i<numClasses; ++i )
      {
      m_SVMClassWeight[i] = 1;
      }
    }
  m_SVMClassWeight[c] = w;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterSVM< TImage, TLabelMap >
::SetSVMClassWeights( VectorDoubleType & w )
{
  unsigned int numClasses = this->m_ObjectIdList.size();

  if( w.size() == numClasses )
    {
    if( m_SVMClassWeight.size() != numClasses )
      {
      m_SVMClassWeight.resize( numClasses );
      }
    for( unsigned int i=0; i<numClasses; ++i )
      {
      m_SVMClassWeight[i] = w[i];
      }
    }
}


template< class TImage, class TLabelMap >
double
PDFSegmenterSVM< TImage, TLabelMap >
::GetSVMClassWeight( unsigned int c ) const
{
  if( c < m_SVMClassWeight.size() )
    {
    return m_SVMClassWeight[ c ];
    }
  return 1;
}

template< class TImage, class TLabelMap >
typename PDFSegmenterSVM< TImage, TLabelMap >::VectorDoubleType &
PDFSegmenterSVM< TImage, TLabelMap >
::GetSVMClassWeights( void )
{
  unsigned int numClasses = this->m_ObjectIdList.size();

  if( m_SVMClassWeight.size() != numClasses )
    {
    m_SVMClassWeight.resize( numClasses );
    for( unsigned int i=0; i<numClasses; ++i )
      {
      m_SVMClassWeight[i] = 1;
      }
    }

  return m_SVMClassWeight;
}

template< class TImage, class TLabelMap >
svm_model *
PDFSegmenterSVM< TImage, TLabelMap >
::GetModel( void )
{
  return m_Model;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterSVM< TImage, TLabelMap >
::SetModel( svm_model * model )
{
  m_Model = model;
}

template< class TImage, class TLabelMap >
svm_parameter *
PDFSegmenterSVM< TImage, TLabelMap >
::GetParameter( void )
{
  return & m_Parameter;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterSVM< TImage, TLabelMap >
::SetParameter( svm_parameter * parameter )
{
  m_Parameter = ( * parameter );
}

template< class TImage, class TLabelMap >
void
PDFSegmenterSVM< TImage, TLabelMap >
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

  m_Parameter.gamma = 1.0 / numFeatures;

  unsigned int sampleSize = 0;
  for( unsigned int c=0; c<numClasses; ++c )
    {
    sampleSize += this->m_InClassList[c].size();
    }
  sampleSize /= m_TrainingDataStride;

  m_Parameter.nr_weight = numClasses;
  m_Parameter.weight_label = ( int * )malloc( numClasses * sizeof( int ) );
  m_Parameter.weight = ( double * )malloc( numClasses * sizeof( double ) );
  if( m_SVMClassWeight.size() != numClasses )
    {
    m_SVMClassWeight.resize( numClasses );
    }
  for( unsigned int c=0; c<numClasses; ++c )
    {
    m_Parameter.weight_label[c] = c;
    m_Parameter.weight[c] = 1.0 - ( 
      ( this->m_InClassList[c].size() / m_TrainingDataStride )
      / ( double )( sampleSize ) );
    m_SVMClassWeight[c] = m_Parameter.weight[c];
    }

  if( m_Problem.l != 0 )
    {
    m_Problem.l = 0;
    if( m_Problem.y != NULL )
      {
      free( m_Problem.y );
      m_Problem.y = NULL;
      }
    if( m_Problem.y != NULL )
      {
      free( m_Problem.y );
      m_Problem.x = NULL;
      }
    }
  m_Problem.l = sampleSize;
  m_Problem.y = ( double * )malloc( sampleSize * sizeof( double ) );
  m_Problem.x = ( struct svm_node ** )malloc( sampleSize *
    sizeof( struct svm_node * ) );

  unsigned int numElements = sampleSize * ( numFeatures + 1 );
  if( m_Space != NULL )
    {
    free( m_Space );
    }
  m_Space = ( struct svm_node * )malloc( numElements *
    sizeof( struct svm_node ) );

  unsigned int elementNum = 0;
  unsigned int sampleNum = 0;
  for( unsigned int c=0; c<numClasses; ++c )
    {
    typename ListSampleType::const_iterator
      inClassListIt( this->m_InClassList[c].begin() );
    typename ListSampleType::const_iterator
      inClassListItEnd( this->m_InClassList[c].end() );
    while( sampleNum < sampleSize && inClassListIt != inClassListItEnd )
      {
      m_Problem.x[ sampleNum ] = & m_Space[ elementNum ];
      m_Problem.y[ sampleNum ] = c;
      for( unsigned int f=0; f<numFeatures; ++f )
        {
        m_Space[ elementNum ].index = f;
        m_Space[ elementNum ].value = ( *inClassListIt )[ f ];
        ++elementNum;
        }
      m_Space[ elementNum ].index = -1;
      m_Space[ elementNum ].value = 0;
      ++elementNum;
      ++sampleNum;
      for( unsigned int s = 0; inClassListIt != inClassListItEnd
        && s < m_TrainingDataStride; ++s )
        {
        ++inClassListIt;
        }
      }
    }

  const char * errorMessage;
  errorMessage = svm_check_parameter( &m_Problem, &m_Parameter );
  if( errorMessage )
    {
    std::cout << "ERROR: LIBSVM = " << errorMessage << std::endl;
    throw;
    }

  m_Model = svm_train( &m_Problem, &m_Parameter );
}

template< class TImage, class TLabelMap >
typename PDFSegmenterSVM< TImage, TLabelMap >::ProbabilityVectorType
PDFSegmenterSVM< TImage, TLabelMap >
::GetProbabilityVector( const FeatureVectorType & fv ) const
{
  unsigned int numClasses = this->m_ObjectIdList.size();
  unsigned int numFeatures = this->GetNumberOfFeatures();

  svm_node * x = ( struct svm_node * )malloc( ( numFeatures+1 ) *
    sizeof( struct svm_node ) );

  for( unsigned int f=0; f<numFeatures; ++f )
    {
    x[ f ].index = f;
    x[ f ].value = fv[ f ];
    }
  x[ numFeatures ].index = -1;
  x[ numFeatures ].value = 0;

  double * probEstimates = ( double * )malloc( numClasses *
    sizeof( double ) );

  svm_predict_probability( m_Model, x, probEstimates );

  ProbabilityVectorType prob( numClasses );
  for( unsigned int c=0; c<numClasses; ++c )
    {
    prob[c] = probEstimates[ c ] * this->GetSVMClassWeight( c );
    }

  free( probEstimates );
  free( x );

  return prob;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterSVM< TImage, TLabelMap >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Training data stride = "
    << m_TrainingDataStride << std::endl;

  if( m_Model == NULL )
    {
    os << indent << "Model = NULL" << std::endl;
    }
  else
    {
    os << indent << "Model = Set" << std::endl;
    }
  os << indent << "Parameter = " << m_Parameter.svm_type << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubePDFSegmenterSVM_hxx )
