/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: main.cxx,v $
  Language:  C++
  Date:      $Date: 2007-07-10 11:35:36 -0400 (Tue, 10 Jul 2007) $
  Version:   $Revision: 48 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkCVTImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"

namespace itk
{

/** Constructor */
template < class TInputImage, class TOutputImage >
CVTImageFilter< TInputImage, TOutputImage >
::CVTImageFilter()
{
  m_Seed = 1234;

  m_NumberOfCentroids = 100;
  m_InitialSamplingMethod = CVT_RANDOM;
  m_NumberOfSamples = 10000;
  m_NumberOfIterations = 100;

  m_BatchSamplingMethod = CVT_RANDOM;
  m_NumberOfIterationsPerBatch = 10;
  m_NumberOfSamplesPerBatch = 5000;
}


/** SetCentroids */
template < class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >
::SetCentroids(const PointArrayType * centroids)
{
  m_InitialSamplingMethod = CVT_USER;
  m_NumberOfCentroids = centroids.size();
  m_Centroids.reserve(m_NumberOfCentroids);
  for(int i = 0; i<m_NumberOfCentroids; i++)
    {
    m_Centroids[i] = centroids[i];
    }
}


/** GenerateInputRequestedRegion */
template < class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
  if ( this->GetInput() )
    {
    typename InputImageType::Pointer inpt = 
          const_cast< TInputImage * >( this->GetInput() );
    inpt->SetRequestedRegionToLargestPossibleRegion();
    }
}


/** EnlargeOutputRequestedRegion */
template < class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >
::EnlargeOutputRequestedRegion(DataObject * output)
{
  Superclass::EnlargeOutputRequestedRegion(output);
  output->SetRequestedRegionToLargestPossibleRegion();
}

/** GenerateData */
template < class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >
::GenerateData()
{
  m_InputImage = this->GetInput();

  RegionType region;
  region = m_InputImage->GetLargestPossibleRegion();

  m_InputImageSize = region.GetSize();

  ImageRegionConstIterator< InputImageType >  inputIt(m_InputImage, region);
  inputIt.GoToBegin();
  m_InputImageMax = inputIt.Get();
  ++inputIt;
  double tf;
  while(!inputIt.IsAtEnd())
    {
    tf = inputIt.Get();
    if( tf > m_InputImageMax )
      {
      m_InputImageMax = tf;
      }
    ++inputIt;
    }
  if(m_InputImageMax == 0)
    {
    m_InputImageMax = 1;
    }

  unsigned int iteration = 0;
  double iterationEnergyDifference = 0.0;
  double iterationEnergy;// = 0.0;

  srand( m_Seed );

  if(m_InitialSamplingMethod != CVT_USER)
    {
    this->ComputeSample( &m_Centroids, m_NumberOfCentroids,
                         m_InitialSamplingMethod );
    }

  m_OutputImage = this->GetOutput(0);

  m_OutputImage->SetRegions(region);
  m_OutputImage->SetOrigin(m_InputImage->GetOrigin());
  m_OutputImage->SetSpacing(m_InputImage->GetSpacing());
  m_OutputImage->Allocate();
  m_OutputImage->FillBuffer( 0 );

  for(int j=0; j<(int)m_NumberOfCentroids; j++)
    {
    std::cout << "Initial Centroid [" << j << "] = " << m_Centroids[j] 
              << std::endl;
    }

  while ( iteration < m_NumberOfIterations )
    {
    iteration = iteration + 1;

    iterationEnergy = this->ComputeIteration(iterationEnergyDifference);

    std::cout << "Iteration = " << iteration 
              << " : E = " << iterationEnergy 
              << " : EDiff = " << iterationEnergyDifference << std::endl;
    ContinuousIndexType indx;
    indx.Fill(0);
    for(int j=0; j<(int)m_NumberOfCentroids; j++)
      {
      for(int i=0; i<ImageDimension; i++)
        {
        indx[i] = indx[i] + m_Centroids[j][i];
        }
      }
    for(int i=0; i<ImageDimension; i++)
      {
      indx[i] = indx[i] / m_NumberOfCentroids;
      }
    std::cout << "   Mean indx = " << indx << std::endl;
    }

  // Generate output image
  IndexType iIndx;
  for(int j=0; j<(int)m_NumberOfCentroids; j++)
    {
    for(int i=0; i<ImageDimension; i++)
      {
      iIndx[i] = (int)(m_Centroids[j][i]);
      if(iIndx[i] < 0)
        {
        iIndx[i] = 0;
        }
      if(iIndx[i] > (int)m_InputImageSize[i]-1)
        {
        iIndx[i] = m_InputImageSize[i]-1;
        }
      }
    m_OutputImage->SetPixel(iIndx, j+1);
    std::cout << " Final Centroid [" << j << "] = " << m_Centroids[j] 
              << std::endl;
    }
  
  typedef DanielssonDistanceMapImageFilter<OutputImageType, OutputImageType>
          DDFilterType;
  typename DDFilterType::Pointer ddFilter = DDFilterType::New();
  ddFilter->SetInput(m_OutputImage);
  ddFilter->SetInputIsBinary(false);
  ddFilter->Update();
  typename OutputImageType::Pointer tmpImage = ddFilter->GetVoronoiMap();

  ImageRegionConstIterator< OutputImageType > tmpImageIt( tmpImage, region );
  ImageRegionIterator< OutputImageType > outputImageIt( m_OutputImage, region );
  tmpImageIt.GoToBegin();
  outputImageIt.GoToBegin();
  while(!tmpImageIt.IsAtEnd())
    {
    outputImageIt.Set(tmpImageIt.Get());
    ++outputImageIt;
    ++tmpImageIt;
    }
}


/** ComputeIteration */
template < class TInputImage, class TOutputImage >
double
CVTImageFilter< TInputImage, TOutputImage >
::ComputeIteration( double & energyDiff)
{
  int i;
  int j;
  int j2;
  double term;

  //  Take each generator as the first sample point for its region.
  //  This can slightly slow the convergence, but it simplifies the
  //  algorithm by guaranteeing that no region is completely missed
  //  by the sampling.
  double energy = 0.0;

  PointArrayType centroids2(m_NumberOfCentroids);
  double * count = new double[m_NumberOfCentroids];
  unsigned int * nearest = new unsigned int[m_NumberOfSamplesPerBatch];
  PointArrayType batch(m_NumberOfSamplesPerBatch);

  for ( j = 0; j < (int)m_NumberOfCentroids; j++ )
    {
    centroids2[j] = m_Centroids[j];
    count[j] = 1;
    }

  //std::cout << " computing iteration..." << std::endl;
  //
  //  Generate the sampling points S.
  //
  int get;
  int have = 0;
  double dist;
  while ( have < (int)m_NumberOfSamples )
    {
    //std::cout << " computing iteration have = " << have << std::endl;
    if( m_NumberOfSamples-have < m_NumberOfSamplesPerBatch )
      {
      get = m_NumberOfSamples - have;
      }
    else
      {
      get = m_NumberOfSamplesPerBatch;
      }

    //std::cout << " computing iteration get = " << get << std::endl;
    ComputeSample( &batch, get, m_BatchSamplingMethod );
    have = have + get;

    ComputeClosest( batch, m_Centroids, nearest );

    for ( j = 0; j < get; j++ )
      {
      j2 = nearest[j];

      dist = 0;
      for ( i = 0; i < ImageDimension; i++ )
        {
        centroids2[j2][i] = centroids2[j2][i] + batch[j][i];
        dist = ( m_Centroids[j2][i] - batch[j][i] )
               * ( m_Centroids[j2][i] - batch[j][i] );
        }
      energy = energy + sqrt(dist);
      count[j2] = count[j2] + 1;
      }
    }

  for ( j = 0; j < (int)m_NumberOfCentroids; j++ )
    {
    for ( i = 0; i < ImageDimension; i++ )
      {
      centroids2[j][i] = centroids2[j][i] / count[j];
      }
    }

  energyDiff = 0.0;
  for ( j = 0; j < (int)m_NumberOfCentroids; j++ )
    {
    term = 0.0;
    for ( i = 0; i < ImageDimension; i++ )
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
template < class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >
::ComputeSample( PointArrayType * sample, unsigned int sampleSize,
                 SamplingMethodEnum samplingMethod )
{
  int i;
  int j;

  if ( sampleSize < 1 )
    {
    std::cout << "\n";
    std::cout << "ComputeSample - Fatal error!\n";
    std::cout << "  sampleSize < 1.\n";
    exit ( 1 );
    }

  (*sample).clear();
  (*sample).reserve(sampleSize);

  //std::cout << "    computing sample" << std::endl;
  //std::cout << "    computing sample size = " << sampleSize << std::endl;
  double p1, u;
  IndexType iIndx;
  ContinuousIndexType indx;
  switch(samplingMethod)
    {
    case CVT_GRID:
      {
      int len = 1;
      for(i=0; i<ImageDimension; i++)
        {
        len = len * m_InputImageSize[i];
        }
      double factor = len / (double)sampleSize;
      factor = pow(factor, 1.0 / ImageDimension);
      int * gridSize = new int[ImageDimension];
      len = 1;
      for(i=0; i<ImageDimension; i++)
        {
        gridSize[i] = (int)(m_InputImageSize[i] / factor);
        len = len * gridSize[i];
        }
      double tmpJ, tmpLen;
      for ( j = 0; j < len; j++ )
        {
        tmpJ = j;
        tmpLen = len;
        for(i=ImageDimension-1; i>=0; i--)
          {
          tmpLen = tmpLen / gridSize[i];
          iIndx[i] = (int)(tmpJ / tmpLen);
          tmpJ = tmpJ - (iIndx[i] * tmpLen);
          iIndx[i] = (int)(iIndx[i] * factor + factor/2);
          }
        (*sample).push_back(iIndx);
        }
      delete [] gridSize;
      for ( j = len; j < (int)sampleSize; j++ )
        {
        for ( i = 0; i < ImageDimension; i++ )
          {
          iIndx[i] = (int)(((double)rand()/(double)RAND_MAX) 
                           * m_InputImageSize[i]);
          if(iIndx[i] < 0)
            {
            iIndx[i] = 0;
            }
          if(iIndx[i] > (int)m_InputImageSize[i]-1)
            {
            iIndx[i] = m_InputImageSize[i]-1;
            }
          }
        (*sample).push_back(iIndx);
        }
      break;
      }
    case CVT_RANDOM:
      {
      for ( j = 0; j < (int)sampleSize; j++ )
        {
        u = 1;
        p1 = 0;
        while(u >= p1)
          {
          for ( i = 0; i < ImageDimension; i++ )
            {
            indx[i] = (((double)rand()/(double)RAND_MAX)*m_InputImageSize[i]);
            iIndx[i] = (int)(indx[i]);
            if(iIndx[i] < 0)
              {
              iIndx[i] = 0;
              }
            if(iIndx[i] > (int)m_InputImageSize[i]-1)
              {
              iIndx[i] = m_InputImageSize[i]-1;
              }
            }
          p1 = m_InputImage->GetPixel(iIndx) / m_InputImageMax;
          u = (double)rand()/(double)RAND_MAX;
          }
        (*sample).push_back(indx);
        }
      break;
      }
    default:
    case CVT_USER:
      {
      std::cout << "\n";
      std::cout << "Sampling method CVT_USER not supported for resampling\n";
      exit ( 1 );
      break;
      }
    }
  //std::cout << "    computing sample done" << std::endl;
}


/** ComputeClosest */
template < class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >
::ComputeClosest( const PointArrayType & sample,
                  const PointArrayType & centroids,
                  unsigned int * nearest )
{
  double distMin;
  double dist;
  int i;
  int jc;
  int js;

  //std::cout << "    computing closest" << std::endl;

  int numberOfSamples = sample.size();

  int numberOfCentroids = centroids.size();

  for ( js = 0; js < numberOfSamples; js++ )
    {
    distMin = 1e20;
    nearest[js] = 0;

    for ( jc = 0; jc < numberOfCentroids; jc++ )
      {
      dist = 0.0;
      for ( i = 0; i < ImageDimension; i++ )
        {
        dist += ( sample[js][i] - centroids[jc][i] ) 
                 * ( sample[js][i] - centroids[jc][i] );
        }

      if ( jc == 0 || dist < distMin )
        {
        distMin = dist;
        nearest[js] = jc;
        }
      }
    }
  // std::cout << "    computing closest done" << std::endl;
}

/** PrintSelf */
template < class TInputImage, class TOutputImage >
void
CVTImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  std::cout << "Input image = " << m_InputImage << std::endl;
  std::cout << "Centroid[0] = " << m_Centroids[0] << std::endl;
  std::cout << "InputImageMax = " << m_InputImageMax << std::endl;
  std::cout << "Seed = " << m_Seed << std::endl;
  std::cout << "NumberOfCentroids = " << m_NumberOfCentroids << std::endl;
  std::cout << "InitialSamplingMethod = " << m_InitialSamplingMethod 
            << std::endl;
  std::cout << "NumberOfSamples = " << m_NumberOfSamples << std::endl;
  std::cout << "NumberOfIterations = " << m_NumberOfIterations << std::endl;
  std::cout << "BatchSamplingMethod = " << m_BatchSamplingMethod << std::endl;
  std::cout << "NumberOfIterationsPerBatch = " << m_NumberOfIterationsPerBatch 
            << std::endl;
  std::cout << "NumberOfSamplesPerBatch = " << m_NumberOfSamplesPerBatch 
            << std::endl;
}

} // namespace itk
