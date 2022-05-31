/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __itktubeTubeExtractorIO_hxx
#define __itktubeTubeExtractorIO_hxx


namespace itk
{

namespace tube
{

template< class TImage >
TubeExtractorIO< TImage >::
TubeExtractorIO( void )
{
  this->Clear();
}

template< class TImage >
TubeExtractorIO< TImage >::
TubeExtractorIO( const char * _headerName )
{
  this->Clear();

  TubeExtractorIO::Read( _headerName );
}

template< class TImage >
TubeExtractorIO< TImage >::
TubeExtractorIO( const typename
  TubeExtractorType::Pointer & _filter )
{
  this->Clear();

  this->InitializeEssential( _filter );
}

template< class TImage >
TubeExtractorIO< TImage >::
~TubeExtractorIO()
{
}

template< class TImage >
void TubeExtractorIO< TImage >::
PrintInfo() const
{
  if( m_TubeExtractor.IsNotNull() )
    {
    std::cout << m_TubeExtractor << std::endl;
    }
  else
    {
    std::cout << "TubeExtractor = NULL" << std::endl;
    }
}

template< class TImage >
void TubeExtractorIO< TImage >::
CopyInfo( const TubeExtractorIO< TImage > & _filterIO )
{
  this->Clear();

  this->InitializeEssential( _filterIO.GetTubeExtractor() );
}

template< class TImage >
void TubeExtractorIO< TImage >::
Clear( void )
{
  m_TubeExtractor = NULL;
}


template< class TImage >
bool TubeExtractorIO< TImage >::
InitializeEssential( const typename
  TubeExtractorType::Pointer & _filter )
{
  m_TubeExtractor = _filter;

  return true;
}

template< class TImage >
void TubeExtractorIO< TImage >::
SetTubeExtractor( TubeExtractorType * _filter )
{
  m_TubeExtractor = _filter;
}

template< class TImage >
const typename TubeExtractor< TImage >::Pointer
TubeExtractorIO< TImage >::
GetTubeExtractor( void ) const
{
  return m_TubeExtractor;
}

template< class TImage >
bool TubeExtractorIO< TImage >::
CanRead( const char * _headerName ) const
{
  MetaTubeExtractor teReader;

  return teReader.CanRead( _headerName );
}

template< class TImage >
bool TubeExtractorIO< TImage >::
Read( const char * _headerName )
{
  if( m_TubeExtractor.IsNull() )
    {
    std::cout <<
      "ERROR: Set a TubeExtractor prior to reading TubeExtractor parameters."
      << std::endl;
    return false;
    }

  typename TubeExtractorType::RidgeExtractorType::Pointer ridgeOp =
    m_TubeExtractor->GetRidgeExtractor();
  typename TubeExtractorType::RadiusExtractorType::Pointer radiusOp =
    m_TubeExtractor->GetRadiusExtractor();

  if( ridgeOp.IsNull() || radiusOp.IsNull() )
    {
    std::cout <<
      "ERROR: Set a tubeExtractor input image prior to reading parameters."
      << std::endl;
    return false;
    }

  MetaTubeExtractor teReader;

  if( !teReader.Read( _headerName ) )
    {
    m_TubeExtractor = NULL;
    return false;
    }

  m_TubeExtractor->SetDataMin( teReader.GetDataMin() );
  m_TubeExtractor->SetDataMax( teReader.GetDataMax() );
  m_TubeExtractor->SetTubeColor( teReader.GetTubeColor() );

  ridgeOp->SetScale( teReader.GetRidgeScale() );
  ridgeOp->SetScaleKernelExtent( teReader.GetRidgeScaleKernelExtent() );
  ridgeOp->SetDynamicScale( teReader.GetRidgeDynamicScale() );
  ridgeOp->SetDynamicStepSize( teReader.GetRidgeDynamicStepSize() );
  ridgeOp->SetStepX( teReader.GetRidgeStepX() );
  ridgeOp->SetMaxTangentChange( teReader.GetRidgeMaxTangentChange() );
  ridgeOp->SetMaxXChange( teReader.GetRidgeMaxXChange() );
  ridgeOp->SetMinRidgeness( teReader.GetRidgeMinRidgeness() );
  ridgeOp->SetMinRidgenessStart( teReader.GetRidgeMinRidgenessStart() );
  ridgeOp->SetMinRoundness( teReader.GetRidgeMinRoundness() );
  ridgeOp->SetMinRoundnessStart( teReader.GetRidgeMinRoundnessStart() );
  ridgeOp->SetMinCurvature( teReader.GetRidgeMinCurvature() );
  ridgeOp->SetMinCurvatureStart( teReader.GetRidgeMinCurvatureStart() );
  ridgeOp->SetMinLevelness( teReader.GetRidgeMinLevelness() );
  ridgeOp->SetMinLevelnessStart( teReader.GetRidgeMinLevelnessStart() );
  ridgeOp->SetMaxRecoveryAttempts( teReader.GetRidgeMaxRecoveryAttempts() );
  ridgeOp->SetDataMin( teReader.GetDataMin() );
  ridgeOp->SetDataMax( teReader.GetDataMax() );

  radiusOp->SetRadiusStart( teReader.GetRadiusStart() );
  radiusOp->SetRadiusMin( teReader.GetRadiusMin() );
  radiusOp->SetRadiusMax( teReader.GetRadiusMax() );
  radiusOp->SetMinMedialness( teReader.GetRadiusMinMedialness() );
  radiusOp->SetMinMedialnessStart( teReader.GetRadiusMinMedialnessStart() );
  radiusOp->SetDataMin( teReader.GetDataMin() );
  radiusOp->SetDataMax( teReader.GetDataMax() );

  return true;
}

template< class TImage >
bool TubeExtractorIO< TImage >::
Write( const char * _headerName )
{
  if( m_TubeExtractor.IsNull() )
    {
    std::cout <<
      "ERROR: Set a tubeExtractor input image prior to writing parameters."
      << std::endl;
    return false;
    }

  MetaTubeExtractor teWriter;

  typename TubeExtractorType::RidgeExtractorType::Pointer ridgeOp =
    m_TubeExtractor->GetRidgeExtractor();
  typename TubeExtractorType::RadiusExtractorType::Pointer radiusOp =
    m_TubeExtractor->GetRadiusExtractor();

  teWriter.SetGeneralProperties( m_TubeExtractor->GetDataMin(),
    m_TubeExtractor->GetDataMax(),
    m_TubeExtractor->GetTubeColor() );

  teWriter.SetRidgeProperties( ridgeOp->GetScale(),
    ridgeOp->GetScaleKernelExtent(),
    ridgeOp->GetDynamicScale(),
    ridgeOp->GetDynamicStepSize(),
    ridgeOp->GetStepX(),
    ridgeOp->GetMaxTangentChange(),
    ridgeOp->GetMaxXChange(),
    ridgeOp->GetMinRidgeness(),
    ridgeOp->GetMinRidgenessStart(),
    ridgeOp->GetMinRoundness(),
    ridgeOp->GetMinRoundnessStart(),
    ridgeOp->GetMinCurvature(),
    ridgeOp->GetMinCurvatureStart(),
    ridgeOp->GetMinLevelness(),
    ridgeOp->GetMinLevelnessStart(),
    ridgeOp->GetMaxRecoveryAttempts() );

  teWriter.SetRadiusProperties( radiusOp->GetRadiusStart(),
    radiusOp->GetRadiusMin(),
    radiusOp->GetRadiusMax(),
    radiusOp->GetMinMedialness(),
    radiusOp->GetMinMedialnessStart() );

  return teWriter.Write( _headerName );
}

} // End namespace tube

} // End namespace itk

#endif // __itktubeTubeExtractorIO_hxx
