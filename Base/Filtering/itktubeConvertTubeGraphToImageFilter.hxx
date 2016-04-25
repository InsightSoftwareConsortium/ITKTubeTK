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
  m_AImage = NULL;
  m_BImage = NULL;
  m_RImage = NULL;
  m_CImage = NULL;
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

  m_AImage = this->GetOutput( 0 );
  m_AImage->SetRegions( m_InputImage->GetLargestPossibleRegion().GetSize() );
  m_AImage->SetSpacing( m_InputImage->GetSpacing() );
  m_AImage->SetOrigin( m_InputImage->GetOrigin() );
  m_AImage->Allocate();
  m_AImage->FillBuffer( 0 );

  m_BImage = OutputImageType::New();
  m_BImage->SetRegions( m_InputImage->GetLargestPossibleRegion().GetSize() );
  m_BImage->SetSpacing( m_InputImage->GetSpacing() );
  m_BImage->SetOrigin( m_InputImage->GetOrigin() );
  m_BImage->Allocate();

  m_RImage = OutputImageType::New();
  m_RImage->SetRegions( m_InputImage->GetLargestPossibleRegion().GetSize() );
  m_RImage->SetSpacing( m_InputImage->GetSpacing() );
  m_RImage->SetOrigin( m_InputImage->GetOrigin() );
  m_RImage->Allocate();

  m_CImage = OutputImageType::New();
  m_CImage->SetRegions( m_InputImage->GetLargestPossibleRegion().GetSize() );
  m_CImage->SetSpacing( m_InputImage->GetSpacing() );
  m_CImage->SetOrigin( m_InputImage->GetOrigin() );
  m_CImage->Allocate();

  itk::ImageRegionConstIterator< InputImageType >
                 itCVT( m_InputImage, m_InputImage->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< OutputImageType >
                 itA( m_AImage, m_AImage->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< OutputImageType >
                 itB( m_BImage, m_BImage->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< OutputImageType >
                 itR( m_RImage, m_RImage->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< OutputImageType >
                 itC( m_CImage, m_CImage->GetLargestPossibleRegion() );
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
