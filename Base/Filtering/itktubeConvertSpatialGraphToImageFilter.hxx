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

#ifndef __itktubeConvertSpatialGraphToImageFilter_hxx
#define __itktubeConvertSpatialGraphToImageFilter_hxx

#include "itktubeConvertSpatialGraphToImageFilter.h"


namespace itk
{

namespace tube
{

/** Constructor */
template< class TInputImage, class TOutputImage >
ConvertSpatialGraphToImageFilter< TInputImage, TOutputImage >::
ConvertSpatialGraphToImageFilter( void )
{
  m_InputImage = NULL;
  m_AdjacencyMatrixImage = NULL;
  m_BranchnessImage = NULL;
  m_RadiusImage = NULL;
  m_CentralityImage = NULL;
}

/** GenerateData */
template< class TInputImage, class TOutputImage >
void
ConvertSpatialGraphToImageFilter< TInputImage, TOutputImage >::
GenerateData( void )
{
  m_InputImage = this->GetInput();
  int numberOfCentroids = m_AdjacencyMatrix.rows();

  m_AdjacencyMatrixImage = this->GetOutput( 0 );
  m_AdjacencyMatrixImage->SetRegions
    ( m_InputImage->GetLargestPossibleRegion().GetSize() );
  m_AdjacencyMatrixImage->SetSpacing( m_InputImage->GetSpacing() );
  m_AdjacencyMatrixImage->SetOrigin( m_InputImage->GetOrigin() );
  m_AdjacencyMatrixImage->Allocate();
  m_AdjacencyMatrixImage->FillBuffer( 0 );

  m_BranchnessImage = OutputImageType::New();
  m_BranchnessImage->SetRegions( m_InputImage->GetLargestPossibleRegion().GetSize() );
  m_BranchnessImage->SetSpacing( m_InputImage->GetSpacing() );
  m_BranchnessImage->SetOrigin( m_InputImage->GetOrigin() );
  m_BranchnessImage->Allocate();

  m_RadiusImage = OutputImageType::New();
  m_RadiusImage->SetRegions( m_InputImage->GetLargestPossibleRegion().GetSize() );
  m_RadiusImage->SetSpacing( m_InputImage->GetSpacing() );
  m_RadiusImage->SetOrigin( m_InputImage->GetOrigin() );
  m_RadiusImage->Allocate();

  m_CentralityImage = OutputImageType::New();
  m_CentralityImage->SetRegions( m_InputImage->GetLargestPossibleRegion().GetSize() );
  m_CentralityImage->SetSpacing( m_InputImage->GetSpacing() );
  m_CentralityImage->SetOrigin( m_InputImage->GetOrigin() );
  m_CentralityImage->Allocate();

  itk::ImageRegionConstIterator< InputImageType >
                 itCVT( m_InputImage, m_InputImage->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< OutputImageType >
                 itA( m_AdjacencyMatrixImage,
                      m_AdjacencyMatrixImage->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< OutputImageType >
                 itB( m_BranchnessImage, m_BranchnessImage->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< OutputImageType >
                 itR( m_RadiusImage, m_RadiusImage->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< OutputImageType >
                 itC( m_CentralityImage, m_CentralityImage->GetLargestPossibleRegion() );
  itCVT.GoToBegin();
  itA.GoToBegin();
  itB.GoToBegin();
  itR.GoToBegin();
  itC.GoToBegin();
  InputPixelType c;

  while( !itCVT.IsAtEnd() )
    {
    c = itCVT.Get()-1;
    itB.Set( m_BranchnessVector[c] );
    itR.Set( m_RadiusVector[c] );
    itC.Set( m_CentralityVector[c] );
    double maxT = 0;
    for( int i = 0; i < numberOfCentroids; i++ )
      {
      if( m_AdjacencyMatrix[c][i] > maxT )
        {
        maxT = m_AdjacencyMatrix[c][i];
        }
      }
    itA.Set( maxT );
    ++itCVT;
    ++itA;
    ++itB;
    ++itR;
    ++itC;
    }
}

template< class TInputImage, class TOutputImage >
void
ConvertSpatialGraphToImageFilter< TInputImage, TOutputImage >::
SetAdjacencyMatrix( vnl_matrix< double > a)
{
  m_AdjacencyMatrix = a;
}

template< class TInputImage, class TOutputImage >
void
ConvertSpatialGraphToImageFilter< TInputImage, TOutputImage >::
SetBranchnessVector( vnl_vector< double > b)
{
  m_BranchnessVector = b;
}

template< class TInputImage, class TOutputImage >
void
ConvertSpatialGraphToImageFilter< TInputImage, TOutputImage >::
SetRadiusVector( vnl_vector< double > r)
{
  m_RadiusVector = r;
}

template< class TInputImage, class TOutputImage >
void
ConvertSpatialGraphToImageFilter< TInputImage, TOutputImage >::
SetCentralityVector( vnl_vector< double > c)
{
  m_CentralityVector = c;
}

/** PrintSelf */
template< class TInputImage, class TOutputImage >
void
ConvertSpatialGraphToImageFilter< TInputImage, TOutputImage >::
PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

} // End namespace tube

} // End namespace itk
#endif // End !defined( __itktubeConvertSpatialGraphToImageFilter_hxx )
