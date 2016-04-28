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

#ifndef __itktubeConvertTubeGraphToImageFilter_hxx
#define __itktubeConvertTubeGraphToImageFilter_hxx

#include "itktubeConvertTubeGraphToImageFilter.h"


namespace itk
{

namespace tube
{

/** Constructor */
template< class TInputImage, class TOutputImage >
ConvertTubeGraphToImageFilter< TInputImage, TOutputImage >::
ConvertTubeGraphToImageFilter( void )
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
ConvertTubeGraphToImageFilter< TInputImage, TOutputImage >::
GenerateData( void )
{
  m_InputImage = this->GetInput();

  double tf;
  int numberOfCentroids;
  int numberOfCentroids2;
  std::ifstream readStream;

  // Read full connectivity information
  std::string matrixFilename = m_InGraphFileName + ".mat";
  readStream.open(matrixFilename.c_str(), std::ios::in);
  readStream >> numberOfCentroids;
  readStream.get();

  vnl_matrix<double> aMat(numberOfCentroids, numberOfCentroids);
  aMat.fill(0);
  vnl_vector<double> bVect(numberOfCentroids);
  bVect.fill(0);
  vnl_vector<double> rVect(numberOfCentroids);
  rVect.fill(0);
  vnl_vector<double> cVect(numberOfCentroids);
  cVect.fill(0);

  for(int i=0; i<numberOfCentroids; i++)
    {
    for(int j=0; j<numberOfCentroids; j++)
      {
      readStream >> tf;
      readStream.get();
      aMat[i][j] = tf;
      }
    }
  readStream.close();

  // Read in BRANCH file
  std::string branchFilename = m_InGraphFileName + ".brc";
  readStream.open(branchFilename.c_str(), std::ios::in);
  readStream >> numberOfCentroids2;
  readStream.get();

  if(numberOfCentroids != numberOfCentroids2)
    {
    return;
    }
  for(int i = 0; i < numberOfCentroids; i++)
    {
    readStream >> tf;
    readStream.get();
    bVect[i] = tf;
    }
  readStream.close();

  // Read in ROOT file
  std::string rootFilename = m_InGraphFileName + ".rot";
  readStream.open(rootFilename.c_str(), std::ios::in);
  readStream >> numberOfCentroids2;
  readStream.get();

  if(numberOfCentroids != numberOfCentroids2)
    {
    return;
    }
  for(int i=0; i<numberOfCentroids; i++)
    {
    readStream >> tf;
    readStream.get();
    rVect[i] = tf;
    }

  // Read in CENTRALITY file
  for(int i=0; i<numberOfCentroids; i++)
    {
    readStream >> tf;
    readStream.get();
    cVect[i] = tf;
    }
  readStream.close();

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
    itB.Set( bVect[c] );
    itR.Set( rVect[c] );
    itC.Set( cVect[c] );
    double maxT = 0;
    for( int i=0; i<numberOfCentroids; i++ )
      {
      if( aMat[c][i] > maxT )
        {
        maxT = aMat[c][i];
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

/** PrintSelf */
template< class TInputImage, class TOutputImage >
void
ConvertTubeGraphToImageFilter< TInputImage, TOutputImage >::
PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeConvertTubeGraphToImageFilter_hxx )
