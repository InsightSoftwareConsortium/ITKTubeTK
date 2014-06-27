/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itktubeRidgeFFTFilter_hxx
#define __itktubeRidgeFFTFilter_hxx

#include "itktubeRidgeFFTFilter.h"

#include "tubeMatrixMath.h"

namespace itk {

namespace tube {
//----------------------------------------------------------------------------
template< typename TInputImage >
RidgeFFTFilter< TInputImage >
::RidgeFFTFilter()
{
  m_Intensity = NULL;
  m_Ridgeness = NULL;
  m_Curvature = NULL;
  m_Levelness = NULL;
  m_Roundness = NULL;

  m_Scale = 1;

  m_DerivativeFilter = DerivativeFilterType::New();
}


template< typename TInputImage >
void
RidgeFFTFilter< TInputImage >
::GenerateData()
{
  std::cout << "Ridge FFT Filter: GenerateData" << std::endl;

  m_DerivativeFilter->SetInput( this->GetInput() );

  typename DerivativeFilterType::OrdersType orders;
  typename DerivativeFilterType::SigmasType sigmas;

  sigmas.Fill( m_Scale );
  m_DerivativeFilter->SetSigmas( sigmas );

  // Intensity
  orders.Fill( 0 );
  m_DerivativeFilter->SetOrders( orders );
  m_DerivativeFilter->Update();
  m_Intensity = m_DerivativeFilter->GetOutput();

  m_Ridgeness = OutputImageType::New();
  m_Ridgeness->CopyInformation( m_Intensity );
  m_Ridgeness->SetRegions( m_Intensity->GetLargestPossibleRegion() );
  m_Ridgeness->Allocate();

  m_Roundness = OutputImageType::New();
  m_Roundness->CopyInformation( m_Intensity );
  m_Roundness->SetRegions( m_Intensity->GetLargestPossibleRegion() );
  m_Roundness->Allocate();

  m_Curvature = OutputImageType::New();
  m_Curvature->CopyInformation( m_Intensity );
  m_Curvature->SetRegions( m_Intensity->GetLargestPossibleRegion() );
  m_Curvature->Allocate();

  m_Levelness = OutputImageType::New();
  m_Levelness->CopyInformation( m_Intensity );
  m_Levelness->SetRegions( m_Intensity->GetLargestPossibleRegion() );
  m_Levelness->Allocate();

  std::vector< typename OutputImageType::Pointer > dx( ImageDimension );
  std::vector< typename OutputImageType::Pointer > ddx( ImageDimension
    * ImageDimension );

  m_DerivativeFilter->GenerateNJet( m_Intensity, dx, ddx );

  ImageRegionIterator< OutputImageType > iterRidge( m_Ridgeness, 
    m_Ridgeness->GetLargestPossibleRegion() );
  ImageRegionIterator< OutputImageType > iterRound( m_Roundness, 
    m_Roundness->GetLargestPossibleRegion() );
  ImageRegionIterator< OutputImageType > iterCurve( m_Curvature, 
    m_Curvature->GetLargestPossibleRegion() );
  ImageRegionIterator< OutputImageType > iterLevel( m_Levelness, 
    m_Levelness->GetLargestPossibleRegion() );

  std::vector< ImageRegionIterator< OutputImageType > > iterDx( 
    ImageDimension );
  std::vector< ImageRegionIterator< OutputImageType > > iterDdx(
    ImageDimension * ImageDimension );
  unsigned int count = 0;
  for( unsigned int i=0; i<ImageDimension; ++i )
    {
    iterDx[i] = ImageRegionIterator< OutputImageType >( dx[i],
      dx[i]->GetLargestPossibleRegion() );
    for( unsigned int j=i; j<ImageDimension; ++j )
      {
      iterDdx[count] = ImageRegionIterator< OutputImageType >( ddx[count],
        ddx[count]->GetLargestPossibleRegion() );
      ++count;
      }
    }

  double ridgeness = 0;
  double roundness = 0;
  double curvature = 0;
  double levelness = 0;
  vnl_matrix<double> H( ImageDimension, ImageDimension);
  vnl_vector<double> D( ImageDimension );
  vnl_matrix<double> HEVect( ImageDimension, ImageDimension);
  vnl_vector<double> HEVal( ImageDimension );
  while( !iterRidge.IsAtEnd() )
    {
    count = 0;
    for( unsigned int i=0; i<ImageDimension; ++i )
      {
      D[i] = iterDx[i].Get();
      ++iterDx[i];
      for( unsigned int j=i; j<ImageDimension; ++j )
        {
        H[i][j] = iterDdx[count].Get();
        H[j][i] = H[i][j];
        ++iterDdx[count];
        ++count;
        }
      }
    ::tube::ComputeRidgeness( H, D, ridgeness, roundness, curvature, levelness,
      HEVect, HEVal );
    iterRidge.Set( ridgeness );
    iterRound.Set( roundness );
    iterCurve.Set( curvature );
    iterLevel.Set( levelness );

    ++iterRidge;
    ++iterRound;
    ++iterCurve;
    ++iterLevel;
    }

  this->SetNthOutput( 0, m_Intensity );
}

template< typename TInputImage >
void
RidgeFFTFilter< TInputImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  if( m_Intensity.IsNotNull() )
    {
    os << indent << "Intensity    : " << m_Intensity << std::endl;
    }
  else
    {
    os << indent << "Intensity    : NULL" << std::endl;
    }

  if( m_Ridgeness.IsNotNull() )
    {
    os << indent << "Ridgeness    : " << m_Ridgeness << std::endl;
    }
  else
    {
    os << indent << "Ridgeness    : NULL" << std::endl;
    }

  if( m_Curvature.IsNotNull() )
    {
    os << indent << "Curvature    : " << m_Curvature << std::endl;
    }
  else
    {
    os << indent << "Curvature    : NULL" << std::endl;
    }

  if( m_Levelness.IsNotNull() )
    {
    os << indent << "Levelness    : " << m_Levelness << std::endl;
    }
  else
    {
    os << indent << "Levelness    : NULL" << std::endl;
    }

  if( m_Roundness.IsNotNull() )
    {
    os << indent << "Roundness    : " << m_Roundness << std::endl;
    }
  else
    {
    os << indent << "Roundness    : NULL" << std::endl;
    }

  os << indent << "Scale             : " << m_Scale << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif
