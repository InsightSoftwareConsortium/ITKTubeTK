/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ITKHeader.h,v $
  Language:  C++
  Date:      $Date: 2007-07-10 11:35:36 -0400 ( Tue, 10 Jul 2007 ) $
  Version:   $Revision: 0 $

  Copyright ( c ) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or https://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itktubeRidgeFFTFilter_hxx
#define __itktubeRidgeFFTFilter_hxx

#include <itkTimeProbesCollectorBase.h>

#include "itktubeFFTGaussianDerivativeIFFTFilter.h"

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
  m_UseIntensityOnly = false;

  m_DerivativeFilter = FFTGaussianDerivativeIFFTFilter< InputImageType,
    OutputImageType >::New();
}


template< typename TInputImage >
void
RidgeFFTFilter< TInputImage >
::GenerateData()
{
  //itk::TimeProbesCollectorBase timeCollector;

  //timeCollector.Start( "RidgeFFT GenerateData" );
  m_DerivativeFilter->SetInput( this->GetInput() );

  typename DerivativeFilterType::OrdersType orders;
  typename DerivativeFilterType::SigmasType sigmas;

  sigmas.Fill( m_Scale );
  m_DerivativeFilter->SetSigmas( sigmas );

  // Intensity
  //timeCollector.Start( "RidgeFFT Intensity" );
  orders.Fill( 0 );
  m_DerivativeFilter->SetOrders( orders );
  m_DerivativeFilter->Update();
  m_Intensity = m_DerivativeFilter->GetOutput();
  //timeCollector.Stop( "RidgeFFT Intensity" );

  if( !m_UseIntensityOnly )
    {
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

    int ddxSize = 0;
    for( unsigned int i=1; i<=ImageDimension; ++i )
      {
      ddxSize += i;
      }
    std::vector< typename OutputImageType::Pointer > ddx( ddxSize );

    //timeCollector.Start( "RidgeFFT GenereateNJet" );
    m_DerivativeFilter->GenerateNJet( m_Intensity, dx, ddx );
    //timeCollector.Stop( "RidgeFFT GenereateNJet" );

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
      ddxSize );

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
    vnl_matrix<double> H( ImageDimension, ImageDimension );
    vnl_vector<double> D( ImageDimension );
    vnl_matrix<double> HEVect( ImageDimension, ImageDimension );
    vnl_vector<double> HEVal( ImageDimension );
    //timeCollector.Start( "RidgeFFT Compute" );
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
      vnl_vector<double> prevTangent;
      ::tube::ComputeRidgeness( H, D, prevTangent, ridgeness, roundness,
        curvature, levelness, HEVect, HEVal );
      iterRidge.Set( ridgeness );
      iterRound.Set( roundness );
      iterCurve.Set( curvature );
      iterLevel.Set( levelness );

      ++iterRidge;
      ++iterRound;
      ++iterCurve;
      ++iterLevel;
      }
    //timeCollector.Stop( "RidgeFFT Compute" );
    }

  this->SetNthOutput( 0, m_Intensity );
  //timeCollector.Stop( "RidgeFFT GenerateData" );
  //timeCollector.Report();
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
  os << indent << "UseIntensityOnly  : " << m_UseIntensityOnly << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif
