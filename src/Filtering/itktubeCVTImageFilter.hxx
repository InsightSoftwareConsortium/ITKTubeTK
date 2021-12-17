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

#ifndef __itktubeCVTImageFilter_hxx
#define __itktubeCVTImageFilter_hxx


#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

namespace itk
{

namespace tube
{

/** Constructor */
template< class TInputImage, class TOutputImage >
CVTImageFilter< TInputImage, TOutputImage >::
CVTImageFilter( void )
{
  m_Seed = -1;
  m_RandomGenerator =
    itk::Statistics::MersenneTwisterRandomVariateGenerator::New();

  m_InputImage = NULL;
  m_InputImageMax = 0;
  m_InputImageSize.Fill( 0 );

  m_OutputImage = NULL;

  m_NumberOfCentroids = 100;
  m_Centroids.clear();

  m_InitialSamplingMethod = CVT_RANDOM;
  m_NumberOfSamples = 10000;
  m_NumberOfIterations = 100;

  m_BatchSamplingMethod = CVT_RANDOM;
  m_NumberOfIterationsPerBatch = 10;
  m_NumberOfSamplesPerBatch = 5000;
}


/** SetCentroids */
template< class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >::
SetCentroids( const PointArrayType & centroids )
{
  m_InitialSamplingMethod = CVT_USER;
  m_NumberOfCentroids = centroids->size();
  m_Centroids = centroids;
}


/** GenerateInputRequestedRegion */
template< class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >::
GenerateInputRequestedRegion( void )
{
  Superclass::GenerateInputRequestedRegion();
  if( this->GetInput() )
    {
    typename InputImageType::Pointer inpt = const_cast< TInputImage * >(
      this->GetInput() );
    inpt->SetRequestedRegionToLargestPossibleRegion();
    }
}


/** EnlargeOutputRequestedRegion */
template< class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >::
EnlargeOutputRequestedRegion( DataObject * output )
{
  Superclass::EnlargeOutputRequestedRegion( output );
  output->SetRequestedRegionToLargestPossibleRegion();
}

/** GenerateData */
template< class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >::
GenerateData( void )
{
  m_InputImage = this->GetInput();

  RegionType region;
  region = m_InputImage->GetLargestPossibleRegion();

  m_InputImageSize = region.GetSize();

  ImageRegionConstIterator< InputImageType > inputIt( m_InputImage,
    region );
  inputIt.GoToBegin();
  m_InputImageMax = inputIt.Get();
  ++inputIt;
  while( !inputIt.IsAtEnd() )
    {
    double tf = inputIt.Get();
    if( tf > m_InputImageMax )
      {
      m_InputImageMax = tf;
      }
    ++inputIt;
    }
  if( m_InputImageMax == 0 )
    {
    m_InputImageMax = 1;
    }

  unsigned int iteration = 0;
  double iterationEnergyDifference = 0.0;

  if( m_Seed != -1 )
    {
    m_RandomGenerator->Initialize( m_Seed );
    }
  else
    {
    m_RandomGenerator->Initialize();
    }

  if( m_InitialSamplingMethod != CVT_USER )
    {
    this->ComputeSample( &m_Centroids, m_NumberOfCentroids,
      m_InitialSamplingMethod );
    }

  m_OutputImage = this->GetOutput( 0 );

  m_OutputImage->SetRegions( region );
  m_OutputImage->SetOrigin( m_InputImage->GetOrigin() );
  m_OutputImage->SetSpacing( m_InputImage->GetSpacing() );
  m_OutputImage->Allocate();
  m_OutputImage->FillBuffer( 0 );

  if( this->GetDebug() )
    {
    for( int j = 0; j < ( int )m_NumberOfCentroids; j++ )
      {
      std::cout << "Initial Centroid [" << j << "] = " << m_Centroids[j]
        << std::endl;
      }
    }

  while( iteration < m_NumberOfIterations )
    {
    iteration = iteration + 1;

    double iterationEnergy = this->ComputeIteration(
      iterationEnergyDifference );

    if( this->GetDebug() )
      {
      std::cout << "Iteration = " << iteration << " : E = "
        << iterationEnergy << " : EDiff = " << iterationEnergyDifference
        << std::endl;
      }
    ContinuousIndexType indx;
    indx.Fill( 0 );
    for( int j = 0; j < ( int )m_NumberOfCentroids; j++ )
      {
      for( unsigned int i=0; i<ImageDimension; i++ )
        {
        indx[i] = indx[i] + m_Centroids[j][i];
        }
      }
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      indx[i] = indx[i] / m_NumberOfCentroids;
      }
    if( this->GetDebug() )
      {
      std::cout << "   Mean indx = " << indx << std::endl;
      }
    }


  // Generate output image
  IndexType iIndx;
  for( int j = 0; j < ( int )m_NumberOfCentroids; j++ )
    {
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      iIndx[i] = ( int )( m_Centroids[j][i] );
      if( iIndx[i] < 0 )
        {
        iIndx[i] = 0;
        }
      if( iIndx[i] > ( int )m_InputImageSize[i]-1 )
        {
        iIndx[i] = m_InputImageSize[i]-1;
        }
      }
    // if you want an output image using the grayscale value at the centroid
    // location, uncomment the following.
    //double tf = m_InputImage->GetPixel( iIndx );
    //m_OutputImage->SetPixel( iIndx, tf ); //j+1 );
    // ...OR...If you want an output imge using the centroid number for each
    // region, uncomment the following.
    m_OutputImage->SetPixel( iIndx, j+1 );
    if( this->GetDebug() )
      {
      std::cout << " Final Centroid [" << j << "] = " << m_Centroids[j]
        << std::endl;
      }
    }

  typedef DanielssonDistanceMapImageFilter<OutputImageType, OutputImageType>
    DDFilterType;
  typename DDFilterType::Pointer ddFilter = DDFilterType::New();
  ddFilter->SetInput( m_OutputImage );
  ddFilter->SetInputIsBinary( false );
  ddFilter->Update();
  typename OutputImageType::Pointer tmpImage = ddFilter->GetVoronoiMap();

  ImageRegionConstIterator< OutputImageType > tmpImageIt( tmpImage,
    region );
  ImageRegionIterator< OutputImageType > outputImageIt( m_OutputImage,
    region );
  tmpImageIt.GoToBegin();
  outputImageIt.GoToBegin();
  while( !tmpImageIt.IsAtEnd() )
    {
    outputImageIt.Set( tmpImageIt.Get() );
    ++outputImageIt;
    ++tmpImageIt;
    }

  this->ComputeAdjacencyMatrix();
}


/** ComputeIteration */
template< class TInputImage, class TOutputImage >
double
CVTImageFilter< TInputImage, TOutputImage >::
ComputeIteration( double & energyDiff )
{
  int i;
  int j;
  int j2;

  //  Take each generator as the first sample point for its region.
  //  This can slightly slow the convergence, but it simplifies the
  //  algorithm by guaranteeing that no region is completely missed
  //  by the sampling.
  double energy = 0.0;

  PointArrayType centroids2( m_NumberOfCentroids );
  double * count = new double[m_NumberOfCentroids];
  unsigned int * nearest = new unsigned int[m_NumberOfSamplesPerBatch];
  PointArrayType batch( m_NumberOfSamplesPerBatch );

  for( j = 0; j < ( int )m_NumberOfCentroids; j++ )
    {
    centroids2[j] = m_Centroids[j];
    count[j] = 1;
    }

  if( this->GetDebug() )
    {
    std::cout << " computing iteration..." << std::endl;
    }
  //
  //  Generate the sampling points S.
  //
  int get;
  int have = 0;
  double dist;
  while( have < ( int )m_NumberOfSamples )
    {
    if( this->GetDebug() )
      {
      std::cout << " computing iteration have = " << have << std::endl;
      }
    if( m_NumberOfSamples-have < m_NumberOfSamplesPerBatch )
      {
      get = m_NumberOfSamples - have;
      }
    else
      {
      get = m_NumberOfSamplesPerBatch;
      }

    if( this->GetDebug() )
      {
      std::cout << " computing iteration get = " << get << std::endl;
      }
    ComputeSample( &batch, get, m_BatchSamplingMethod );
    have = have + get;

    ComputeClosest( batch, m_Centroids, nearest );

    for( j = 0; j < get; j++ )
      {
      j2 = nearest[j];

      dist = 0;
      for( i = 0; i < ( int )ImageDimension; i++ )
        {
        centroids2[j2][i] = centroids2[j2][i] + batch[j][i];
        dist = ( m_Centroids[j2][i] - batch[j][i] )
          * ( m_Centroids[j2][i] - batch[j][i] );
        }
      energy = energy + std::sqrt( dist );
      count[j2] = count[j2] + 1;
      }
    }

  for( j = 0; j < ( int )m_NumberOfCentroids; j++ )
    {
    for( i = 0; i < ( int )ImageDimension; i++ )
      {
      centroids2[j][i] = centroids2[j][i] / count[j];
      }
    }

  energyDiff = 0.0;
  for( j = 0; j < ( int )m_NumberOfCentroids; j++ )
    {
    double term = 0.0;
    for( i = 0; i < ( int )ImageDimension; i++ )
      {
      term += ( centroids2[j][i] - m_Centroids[j][i] )
        * ( centroids2[j][i] - m_Centroids[j][i] );
      m_Centroids[j][i] = centroids2[j][i];
      }
    energyDiff += sqrt ( term );
    }

  energy = energy / m_NumberOfSamples;

  delete [] count;
  delete [] nearest;

  return energy;
}


/** ComputeSample */
template< class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >::
ComputeSample( PointArrayType * sample, unsigned int sampleSize,
  SamplingMethodEnum samplingMethod )
{
  if( sampleSize < 1 )
    {
    if( this->GetDebug() )
      {
      std::cout << "ComputeSample: Fatal error! SampleSize < 1."
        << std::endl;
      }
    throw( "Sample size less than 1.  Cannot compute CVT" );
    }

  sample->clear();
  sample->reserve( sampleSize );

  if( this->GetDebug() )
    {
    std::cout << "    computing sample" << std::endl;
    std::cout << "    computing sample size = " << sampleSize << std::endl;
    }
  IndexType iIndx;
  iIndx.Fill( 0 );
  ContinuousIndexType indx;
  switch( samplingMethod )
    {
    case CVT_GRID:
      {
      unsigned int len = 1;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        len = len * m_InputImageSize[i];
        }
      double factor = len / ( double )sampleSize;
      factor = std::pow( factor, 1.0 / ImageDimension );
      int * gridSize = new int[ImageDimension];
      len = 1;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        gridSize[i] = ( int )( m_InputImageSize[i] / factor );
        len = len * gridSize[i];
        }
      for( unsigned int j = 0; j < len; j++ )
        {
        double tmpJ = j;
        double tmpLen = len;
        for( int i = ( int )( ImageDimension )-1; i >= 0; --i )
          {
          tmpLen = tmpLen / gridSize[i];
          iIndx[i] = ( int )( tmpJ / tmpLen );
          tmpJ = tmpJ - ( iIndx[i] * tmpLen );
          iIndx[i] = ( int )( iIndx[i] * factor + factor/2 );
          }
        sample->push_back( iIndx );
        }
      delete [] gridSize;
      for( unsigned int j = len; j < sampleSize; j++ )
        {
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          iIndx[i] = ( int )( m_RandomGenerator->GetUniformVariate( 0, 1 )
            * m_InputImageSize[i]-1 );
          }
        ( *sample ).push_back( iIndx );
        }
      break;
      }
    case CVT_RANDOM:
      {
      for( unsigned int j = 0; j < sampleSize; j++ )
        {
        double u = 1;
        double p1 = 0;
        while( u >= p1 )
          {
          for( unsigned int i = 0; i < ImageDimension; i++ )
            {
            indx[i] = ( int )( m_RandomGenerator->GetUniformVariate( 0, 1 )
              * m_InputImageSize[i]-1 );
            iIndx[i] = ( int )( indx[i] );
            }
          p1 = m_InputImage->GetPixel( iIndx ) / m_InputImageMax;
          u = ( double )m_RandomGenerator->GetUniformVariate( 0, 1 );
          }
        sample->push_back( indx );
        }
      break;
      }
    default:
    case CVT_USER:
      {
      if( this->GetDebug() )
        {
        std::cout << "Sampling method CVT_USER not supported for resampling"
          << std::endl;
        }
      throw( "Sampling method CVT_USER not supported for resmpling." );
      }
    }
}


/** ComputeClosest */
template< class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >::
ComputeClosest( const PointArrayType & sample,
  const PointArrayType & centroids, unsigned int * nearest )
{
  if( this->GetDebug() )
    {
    std::cout << "    computing closest" << std::endl;
    }

  unsigned int numberOfSamples = sample.size();

  unsigned int numberOfCentroids = centroids.size();

  for( unsigned int js = 0; js < numberOfSamples; js++ )
    {
    double distMin = 1e20;
    nearest[js] = 0;

    for( unsigned int jc = 0; jc < numberOfCentroids; jc++ )
      {
      double dist = 0.0;
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        dist += ( sample[js][i] - centroids[jc][i] )
          * ( sample[js][i] - centroids[jc][i] );
        }

      if( jc == 0 || dist < distMin )
        {
        distMin = dist;
        nearest[js] = jc;
        }
      }
    }
  if( this->GetDebug() )
    {
    std::cout << "    computing closest done" << std::endl;
    }
}

template< class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >::
ComputeAdjacencyMatrix( void )
{
  m_AdjacencyMatrix.SetSize( m_NumberOfCentroids, m_NumberOfCentroids );
  m_AdjacencyMatrix.Fill( 0 );

  typename OutputImageType::SizeType size =
    m_OutputImage->GetLargestPossibleRegion().GetSize();
  Index<ImageDimension> indx;
  Index<ImageDimension> indx2;
  Index<ImageDimension> indx3;
  indx.Fill( 0 );
  bool done = false;
  int n;
  while( !done )
    {
    int c = ( int )( m_OutputImage->GetPixel( indx )-1 );
    indx2.Fill( 0 );
    indx2[0] = 1;
    bool done2 = false;
    while( !done2 )
      {
      bool invalid = false;
      for( unsigned int d=0; d<ImageDimension; d++ )
        {
        indx3[d] = indx[d] + indx2[d];
        if( indx3[d] >= ( int )size[d] )
          {
          invalid = true;
          break;
          }
        }
      if( !invalid )
        {
        n = ( int )( m_OutputImage->GetPixel( indx3 )-1 );
        if( c != n )
          {
          m_AdjacencyMatrix( c, n ) = 1;
          m_AdjacencyMatrix( n, c ) = 1;
          }
        }
      int i=0;
      indx2[i]++;
      while( !done2 && indx2[i]>=2 )
        {
        indx2[i] = 0;
        i++;
        if( i > ( int )( ImageDimension )-1 )
          {
          done2 = true;
          }
        else
          {
          indx2[i]++;
          }
        }
      }
    int i = 0;
    indx[i]++;
    while( !done && indx[i]>=( int )size[i] )
      {
      indx[i] = 0;
      i++;
      if( i>( int )( ImageDimension )-1 )
        {
        done = true;
        }
      else
        {
        if( i == ( int )( ImageDimension )-1 )
          {
          std::cout << "Computing adjacency of slice : "
                    << indx[ImageDimension-1] << std::endl;
          }
        indx[i]++;
        }
      }
    }
}

/** PrintSelf */
template< class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >::
PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_InputImage.IsNotNull() )
    {
    std::cout << "Input image = " << m_InputImage << std::endl;
    }
  else
    {
    std::cout << "Input image = NULL" << std::endl;
    }
  std::cout << "InputImageMax = " << m_InputImageMax << std::endl;
  std::cout << "InputImageSize = " << m_InputImageSize << std::endl;

  if( m_OutputImage.IsNotNull() )
    {
    std::cout << "Output image = " << m_OutputImage << std::endl;
    }
  else
    {
    std::cout << "Output image = NULL" << std::endl;
    }

  std::cout << "NumberOfCentroids = " << m_NumberOfCentroids << std::endl;
  if( m_Centroids.size() > 0 )
    {
    std::cout << "Centroid[0] = " << m_Centroids[0] << std::endl;
    }
  else
    {
    std::cout << "Centroid = NULL" << std::endl;
    }

  std::cout << "Seed = " << m_Seed << std::endl;
  std::cout << "InitialSamplingMethod = " << m_InitialSamplingMethod
    << std::endl;
  std::cout << "NumberOfSamples = " << m_NumberOfSamples << std::endl;
  std::cout << "NumberOfIterations = " << m_NumberOfIterations << std::endl;
  std::cout << "BatchSamplingMethod = " << m_BatchSamplingMethod
    << std::endl;
  std::cout << "NumberOfIterationsPerBatch = "
    << m_NumberOfIterationsPerBatch << std::endl;
  std::cout << "NumberOfSamplesPerBatch = " << m_NumberOfSamplesPerBatch
    << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeCVTImageFilter_hxx )
