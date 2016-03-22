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
  m_Model = NULL;

  m_Parameter.svm_type = C_SVC;
  m_Parameter.kernel_type = RBF;
  m_Parameter.degree = 3;
  m_Parameter.gamma = 0; // 1/num_features
  m_Parameter.coef0 = 0;
  m_Parameter.nu = 0.5;
  m_Parameter.cache_size = 100;
  m_Parameter.C = 1;
  m_Parameter.eps = 1e-3;
  m_Parameter.p = 0.1;
  m_Parameter.shrinking = 1;
  m_Parameter.probability = 1;
  m_Parameter.nr_weight = 0;
  m_Parameter.weight_label = NULL;
  m_Parameter.weight = NULL;
}

template< class TImage, class TLabelMap >
PDFSegmenterSVM< TImage, TLabelMap >
::~PDFSegmenterSVM( void )
{
  if( m_Model != NULL )
    {
    svm_free_and_destroy_model( & m_Model );
    }
  svm_destroy_param( & m_Parameter );
}

template< class TImage, class TLabelMap >
svm_model *
PDFSegmenterSVM< TImage, TLabelMap >
::GetModel( void )
{
  return m_Model;
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

  m_Parameter.nr_weight = numClasses;
  m_Parameter.weight_label = (int *)malloc( numClasses * sizeof( int ) );
  m_Parameter.weight = (double *)malloc( numClasses * sizeof( double ) );
  for( unsigned int c=0; c<numClasses; ++c )
    {
    m_Parameter.weight_label[c] = c;
    m_Parameter.weight[c] = 1.0 - ( this->m_InClassList[c].size() /
      (double)( sampleSize ) );
    }

  svm_problem prob;
  prob.l = sampleSize;
  prob.y = (double *)malloc( sampleSize * sizeof( double ) );
  prob.x = (struct svm_node **)malloc( sampleSize *
    sizeof( struct svm_node *) );

  unsigned int numElements = sampleSize * ( numFeatures + 1 );
  svm_node *x_space;
  x_space = (struct svm_node *)malloc( numElements *
    sizeof( struct svm_node ) );

  unsigned int elementNum = 0;
  unsigned int sampleNum = 0;
  for( unsigned int c=0; c<numClasses; ++c )
    {
    typename ListSampleType::const_iterator
      inClassListIt( this->m_InClassList[c].begin() );
    typename ListSampleType::const_iterator
      inClassListItEnd( this->m_InClassList[c].end() );
    while( inClassListIt != inClassListItEnd )
      {
      prob.x[ sampleNum ] = & x_space[ elementNum ];
      prob.y[ sampleNum ] = c;
      for( unsigned int f=0; f<numFeatures; ++f )
        {
        x_space[ elementNum ].index = f;
        x_space[ elementNum ].value = (*inClassListIt)[ f ];
        ++elementNum;
        }
      x_space[ elementNum ].index = -1;
      x_space[ elementNum ].value = 0;
      ++elementNum;
      ++sampleNum;
      ++inClassListIt;
      }
    }
  const char * errorMessage;
  errorMessage = svm_check_parameter( &prob, &m_Parameter );
  if( errorMessage )
    {
    std::cout << "ERROR: LIBSVM = " << errorMessage << std::endl;
    throw;
    }

  m_Model = svm_train( &prob, &m_Parameter );

  //free( prob.y );
  //free( prob.x );
  //free( x_space );
}

template< class TImage, class TLabelMap >
typename PDFSegmenterSVM< TImage, TLabelMap >::ProbabilityPixelType
PDFSegmenterSVM< TImage, TLabelMap >
::GetClassProbability( unsigned int classNum, const FeatureVectorType & fv)
  const
{
  unsigned int numClasses = this->m_ObjectIdList.size();
  unsigned int numFeatures = this->GetNumberOfFeatures();

  svm_node * x = (struct svm_node *)malloc( ( numFeatures+1 ) *
    sizeof( struct svm_node ) );

  for( unsigned int f=0; f<numFeatures; ++f )
    {
    x[ f ].index = f;
    x[ f ].value = fv[ f ];
    }
  x[ numFeatures ].index = -1;
  x[ numFeatures ].value = 0;

  double * probEstimates = (double *)malloc( numClasses *
    sizeof( double ) );

  double classEstimate = svm_predict_probability( m_Model, x,
    probEstimates );
  double classProb = probEstimates[ classNum ];

  free( probEstimates );
  free( x );

  return classProb;
}

template< class TImage, class TLabelMap >
void
PDFSegmenterSVM< TImage, TLabelMap >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  //os << indent << "Model = " << m_Model << std::endl;
  //os << indent << "Parameter = " << m_Parameter << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubePDFSegmenterSVM_hxx)
