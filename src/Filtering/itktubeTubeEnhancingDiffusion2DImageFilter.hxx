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

#ifndef __itktubeTubeEnhancingDiffusion2DImageFilter_hxx
#define __itktubeTubeEnhancingDiffusion2DImageFilter_hxx

#include "itktubeTubeEnhancingDiffusion2DImageFilter.h"

#include <itkCastImageFilter.h>
#include <itkConstShapedNeighborhoodIterator.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkMinimumMaximumImageFilter.h>
#include <itkNeighborhoodAlgorithm.h>
#include <itkNumericTraits.h>
#include <itkProgressAccumulator.h>
#include <itkZeroFluxNeumannBoundaryCondition.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <iostream>

namespace itk
{

namespace tube
{

template< class TPixel, unsigned int TDimension >
TubeEnhancingDiffusion2DImageFilter<TPixel, TDimension>
::TubeEnhancingDiffusion2DImageFilter( void )
  : m_TimeStep( 0.2 ),
    m_Iterations( 200 ),
    m_RecalculateTubeness( 100 ),
    m_Beta( 0.5 ),
    m_Gamma( 5.0 ),
    m_Epsilon( 0.001 ),
    m_Omega( 25.0 ),
    m_Sensitivity( 20.0 ),
    m_DarkObjectLightBackground( false ),
    m_Verbose( false ),
    m_CurrentIteration( 0 )
{
  this->SetNumberOfRequiredInputs( 1 );

  m_Scales.resize( 2 );
  m_Scales[0] = 6;
  m_Scales[1] = 8;
}


template< class TPixel, unsigned int TDimension >
void
TubeEnhancingDiffusion2DImageFilter<TPixel, TDimension>
::PrintSelf( std::ostream &os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "TimeStep                  : " << m_TimeStep << std::endl;
  os << indent << "Iterations                : " << m_Iterations
    << std::endl;
  os << indent << "RecalculateTubeness       : " << m_RecalculateTubeness
     << std::endl;
  os << indent << "Scales                    : ";
  for( unsigned int i=0; i<m_Scales.size(); ++i )
    {
    os << m_Scales[i] << " ";
    }
  os << std::endl;
  os << indent << "Epsilon                   : " << m_Epsilon << std::endl;
  os << indent << "Omega                     : " << m_Omega << std::endl;
  os << indent << "Sensitivity               : " << m_Sensitivity
    << std::endl;
  os << indent << "DarkObjectLightBackground : "
    << m_DarkObjectLightBackground << std::endl;
  os << indent << "Beta                      : " << m_Beta << std::endl;
  os << indent << "Gamma                     : " << m_Gamma << std::endl;
  os << indent << "Verbose                   : " << m_Verbose << std::endl;
}


template< class TPixel, unsigned int TDimension >
void
TubeEnhancingDiffusion2DImageFilter<TPixel, TDimension>
::VED2DSingleIteration( typename PrecisionImageType::Pointer ci )
{
  bool rec( false );
  if( ( m_CurrentIteration == 1 ) ||
       ( m_RecalculateTubeness == 0 ) ||
       ( m_CurrentIteration % m_RecalculateTubeness == 0 ) )
    {
    rec = true;
    if( m_Verbose )
      {
      std::cout << "v ";
      std::cout.flush();
      }
    MaxTubeResponse( ci );
    DiffusionTensor();
    }
  if( m_Verbose )
    {
    if( !rec )
      {
      std::cout << ". ";
      std::cout.flush();
      }
    }

  // calculate d = nonlineardiffusion( ci )
  // using 3x3x3 stencil, afterwards copy
  // result from d back to ci
  typename PrecisionImageType::Pointer d = PrecisionImageType::New();
  d->SetOrigin( ci->GetOrigin() );
  d->SetSpacing( ci->GetSpacing() );
  d->SetDirection( ci->GetDirection() );
  d->SetRegions( ci->GetLargestPossibleRegion() );
  d->Allocate();
  d->FillBuffer( NumericTraits<Precision>::Zero );

  // shapedneighborhood iter, zeroflux boundary condition
  // division into faces and inner region
  typedef ZeroFluxNeumannBoundaryCondition<PrecisionImageType>     BT;
  typedef ConstShapedNeighborhoodIterator<PrecisionImageType, BT>  NT;
  typedef typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<
    PrecisionImageType>                                            FT;

  BT                      b;
  typename NT::RadiusType r;
  r.Fill( 1 );

  // offsets
  const typename NT::OffsetType oxp = {{1, 0}};
  const typename NT::OffsetType oxm = {{-1, 0}};
  const typename NT::OffsetType oyp = {{0, 1}};
  const typename NT::OffsetType oym = {{0, -1}};

  const typename NT::OffsetType oxpyp = {{1, 1}};
  const typename NT::OffsetType oxmym = {{-1, -1}};
  const typename NT::OffsetType oxpym = {{1, -1}};
  const typename NT::OffsetType oxmyp = {{-1, 1}};

  // fixed weights ( timers )
  const typename PrecisionImageType::SpacingType ispacing =
    ci->GetSpacing();
  const Precision rxx = m_TimeStep / ( 2.0 * ispacing[0] * ispacing[0] );
  const Precision ryy = m_TimeStep / ( 2.0 * ispacing[1] * ispacing[1] );
  const Precision rxy = m_TimeStep / ( 4.0 * ispacing[0] * ispacing[1] );

  // faces
  FT                        fc;
  typename FT::FaceListType fci = fc( ci, d->GetLargestPossibleRegion(),
    r );
  typename FT::FaceListType fxx = fc( m_Dxx, d->GetLargestPossibleRegion(),
    r );
  typename FT::FaceListType fxy = fc( m_Dxy, d->GetLargestPossibleRegion(),
    r );
  typename FT::FaceListType fyy = fc( m_Dyy, d->GetLargestPossibleRegion(),
    r );

  typename FT::FaceListType::iterator fitci, fitxx, fitxy, fityy;

  for( fitci = fci.begin(), fitxx = fxx.begin(), fitxy = fxy.begin(),
    fityy = fyy.begin(); fitci != fci.end();
    ++fitci, ++fitxx, ++fitxy, ++fityy )
    {
    // output iter
    ImageRegionIterator<PrecisionImageType> dit( d, *fitci );

    // input iters
    NT itci( r, ci, *fitci );
    NT itxx( r, m_Dxx, *fitxx );
    NT itxy( r, m_Dxy, *fitxy );
    NT ityy( r, m_Dyy, *fityy );

    itci.OverrideBoundaryCondition( &b );
    itxx.OverrideBoundaryCondition( &b );
    itxy.OverrideBoundaryCondition( &b );
    ityy.OverrideBoundaryCondition( &b );

    // setting active offsets ( yeah there must
    // be some smarter way of doing this )
    itci.ClearActiveList();
    itxx.ClearActiveList();
    itxy.ClearActiveList();
    ityy.ClearActiveList();

    itci.ActivateOffset( oxp );
    itxx.ActivateOffset( oxp );
    itxy.ActivateOffset( oxp );
    ityy.ActivateOffset( oxp );

    itci.ActivateOffset( oxm );
    itxx.ActivateOffset( oxm );
    itxy.ActivateOffset( oxm );
    ityy.ActivateOffset( oxm );

    itci.ActivateOffset( oyp );
    itxx.ActivateOffset( oyp );
    itxy.ActivateOffset( oyp );
    ityy.ActivateOffset( oyp );

    itci.ActivateOffset( oym );
    itxx.ActivateOffset( oym );
    itxy.ActivateOffset( oym );
    ityy.ActivateOffset( oym );

    itci.ActivateOffset( oxpyp );
    itxx.ActivateOffset( oxpyp );
    itxy.ActivateOffset( oxpyp );
    ityy.ActivateOffset( oxpyp );

    itci.ActivateOffset( oxmym );
    itxx.ActivateOffset( oxmym );
    itxy.ActivateOffset( oxmym );
    ityy.ActivateOffset( oxmym );

    itci.ActivateOffset( oxpym );
    itxx.ActivateOffset( oxpym );
    itxy.ActivateOffset( oxpym );
    ityy.ActivateOffset( oxpym );

    itci.ActivateOffset( oxmyp );
    itxx.ActivateOffset( oxmyp );
    itxy.ActivateOffset( oxmyp );
    ityy.ActivateOffset( oxmyp );

    // run for each face diffusion
    for( itci.GoToBegin(), dit.GoToBegin(), itxx.GoToBegin(),
         itxy.GoToBegin(), ityy.GoToBegin();
         !itci.IsAtEnd();
         ++itci, ++dit, ++itxx, ++itxy, ++ityy )
      {
      // weights
      const Precision xp = itxx.GetPixel( oxp ) + itxx.GetCenterPixel();
      const Precision xm = itxx.GetPixel( oxm ) + itxx.GetCenterPixel();
      const Precision yp = ityy.GetPixel( oyp ) + ityy.GetCenterPixel();
      const Precision ym = ityy.GetPixel( oym ) + ityy.GetCenterPixel();

      const Precision xpyp = itxy.GetPixel( oxpyp ) + itxy.GetCenterPixel();
      const Precision xmym = itxy.GetPixel( oxmym ) + itxy.GetCenterPixel();
      const Precision xpym = - itxy.GetPixel( oxpym )
        - itxy.GetCenterPixel();
      const Precision xmyp = - itxy.GetPixel( oxmyp )
        - itxy.GetCenterPixel();

      // evolution
      const Precision cv = itci.GetCenterPixel();
      dit.Value() = cv +
            + rxx * ( xp * ( itci.GetPixel( oxp ) - cv )
                    + xm * ( itci.GetPixel( oxm ) - cv ) )
            + ryy * ( yp * ( itci.GetPixel( oyp ) - cv )
                    + ym * ( itci.GetPixel( oym ) - cv ) )
            + rxy * ( xpyp * ( itci.GetPixel( oxpyp ) - cv )
                    + xmym * ( itci.GetPixel( oxmym ) - cv )
                    + xpym * ( itci.GetPixel( oxpym ) - cv )
                    + xmyp * ( itci.GetPixel( oxmyp ) - cv ) );

      }
    }

  // copying
  ImageRegionConstIterator<PrecisionImageType>
    iti ( d, d->GetLargestPossibleRegion() );
  ImageRegionIterator<PrecisionImageType>
    ito ( ci, ci->GetLargestPossibleRegion() );
  for( iti.GoToBegin(), ito.GoToBegin(); !iti.IsAtEnd(); ++iti, ++ito )
    {
    ito.Value() = iti.Value();
    }

  return;
}


template< class TPixel, unsigned int TDimension >
void
TubeEnhancingDiffusion2DImageFilter<TPixel, TDimension>
::MaxTubeResponse( const typename PrecisionImageType::Pointer im )
{

  // alloc memory for Hessian/tensor
  m_Dxx = PrecisionImageType::New();
  m_Dxx->SetOrigin( im->GetOrigin() );
  m_Dxx->SetSpacing( im->GetSpacing() );
  m_Dxx->SetDirection( im->GetDirection() );
  m_Dxx->SetRegions( im->GetLargestPossibleRegion() );
  m_Dxx->Allocate();
  m_Dxx->FillBuffer( NumericTraits<Precision>::One );

  m_Dxy = PrecisionImageType::New();
  m_Dxy->SetOrigin( im->GetOrigin() );
  m_Dxy->SetSpacing( im->GetSpacing() );
  m_Dxy->SetDirection( im->GetDirection() );
  m_Dxy->SetRegions( im->GetLargestPossibleRegion() );
  m_Dxy->Allocate();
  m_Dxy->FillBuffer( NumericTraits<Precision>::Zero );

  m_Dyy = PrecisionImageType::New();
  m_Dyy->SetOrigin( im->GetOrigin() );
  m_Dyy->SetSpacing( im->GetSpacing() );
  m_Dyy->SetDirection( im->GetDirection() );
  m_Dyy->SetRegions( im->GetLargestPossibleRegion() );
  m_Dyy->Allocate();
  m_Dyy->FillBuffer( NumericTraits<Precision>::One );

  // create temp vesselness image to store maxvessel
  typename PrecisionImageType::Pointer vi = PrecisionImageType::New();

  vi->SetOrigin( im->GetOrigin() );
  vi->SetSpacing( im->GetSpacing() );
  vi->SetDirection( im->GetDirection() );
  vi->SetRegions( im->GetLargestPossibleRegion() );
  vi->Allocate();
  vi->FillBuffer( NumericTraits<Precision>::Zero );

  for( unsigned int i = 0; i < m_Scales.size(); ++i )
    {
    typedef HessianRecursiveGaussianImageFilter<
      PrecisionImageType > HessianType;
    typename HessianType::Pointer hessian = HessianType::New();
    hessian->SetInput( im );
    hessian->SetNormalizeAcrossScale( true );
    hessian->SetSigma( m_Scales[i] );
    hessian->Update();

    ImageRegionIterator<PrecisionImageType>
      itxx ( m_Dxx, m_Dxx->GetLargestPossibleRegion() );
    ImageRegionIterator<PrecisionImageType>
      itxy ( m_Dxy, m_Dxy->GetLargestPossibleRegion() );
    ImageRegionIterator<PrecisionImageType>
      ityy ( m_Dyy, m_Dyy->GetLargestPossibleRegion() );
    ImageRegionIterator<PrecisionImageType>
      vit( vi, vi->GetLargestPossibleRegion() );

    ImageRegionConstIterator<typename HessianType::OutputImageType> hit(
      hessian->GetOutput(),
      hessian->GetOutput()->GetLargestPossibleRegion() );

    for( itxx.GoToBegin(), itxy.GoToBegin(), ityy.GoToBegin(),
         vit.GoToBegin(), hit.GoToBegin();
         !vit.IsAtEnd();
         ++itxx, ++itxy, ++ityy, ++hit, ++vit )
      {
      vnl_matrix<Precision> H( 2, 2 );

      H( 0, 0 ) = hit.Value()( 0, 0 );
      H( 0, 1 ) = H( 1, 0 ) = hit.Value()( 0, 1 );
      H( 1, 1 ) = hit.Value()( 1, 1 );

      vnl_symmetric_eigensystem<Precision> ES( H );
      vnl_vector<Precision> ev( 2 );

      ev[0] = ES.get_eigenvalue( 0 );
      ev[1] = ES.get_eigenvalue( 1 );

      if( std::fabs( ev[0] ) > std::fabs( ev[1] ) )
        {
        std::swap( ev[0], ev[1] );
        }

      const Precision vesselness = TubenessFunction2D( ev[0], ev[1] );

      if( vesselness > 0 && vesselness > vit.Value() )
        {
        vit.Value() = vesselness;

        itxx.Value() = hit.Value()( 0, 0 );
        itxy.Value() = hit.Value()( 0, 1 );
        ityy.Value() = hit.Value()( 1, 1 );
        }
      }
    }

  return;
}


template< class TPixel, unsigned int TDimension >
typename TubeEnhancingDiffusion2DImageFilter<TPixel, TDimension>::Precision
TubeEnhancingDiffusion2DImageFilter<TPixel, TDimension>
::TubenessFunction2D( const Precision l1, const Precision l2 )
{
  if( l2 < 0 )
    {
    return 0;
    }

  Precision vesselness;

  //const Precision smoothC=1E-5;

  const Precision vb2= 2.0*m_Beta*m_Beta;
  const Precision vc2= 2.0*m_Gamma*m_Gamma;

  const Precision   Rb2 = ( l1 * l1 ) / ( l2 * l2 ); // btwn 0 and 1
  const Precision   S2 =  ( l1 * l1 ) + ( l2 *l2 );
  //const Precision   T = std::exp( -( 2*smoothC*smoothC )
  //  / ( std::fabs( l1 )*l2*l2 ) );

  vesselness = std::exp( -Rb2/vb2 ) * ( 1.0 - std::exp( -S2/vc2 ) );

  return vesselness;
}


template< class TPixel, unsigned int TDimension >
void
TubeEnhancingDiffusion2DImageFilter<TPixel, TDimension>
::DiffusionTensor( void )
{
  ImageRegionIterator<PrecisionImageType>
    itxx( m_Dxx, m_Dxx->GetLargestPossibleRegion() );
  ImageRegionIterator<PrecisionImageType>
    itxy( m_Dxy, m_Dxy->GetLargestPossibleRegion() );
  ImageRegionIterator<PrecisionImageType>
    ityy( m_Dyy, m_Dyy->GetLargestPossibleRegion() );

  for( itxx.GoToBegin(), itxy.GoToBegin(), ityy.GoToBegin();
       !itxx.IsAtEnd();
       ++itxx, ++itxy, ++ityy )
    {
    vnl_matrix<Precision> H( 2, 2 );
    H( 0, 0 ) = itxx.Value();
    H( 0, 1 ) = H( 1, 0 ) = itxy.Value();
    H( 1, 1 ) = ityy.Value();

    vnl_symmetric_eigensystem<Precision> ES( H );

    vnl_matrix<Precision> EV( 2, 2 );
    EV.set_column( 0, ES.get_eigenvector( 0 ) );
    EV.set_column( 1, ES.get_eigenvector( 1 ) );

    vnl_vector<Precision> ev( 2 );
    ev[0] = ES.get_eigenvalue( 0 );
    ev[1] = ES.get_eigenvalue( 1 );

    if( std::fabs( ev[0] ) > std::fabs( ev[1] ) )
      {
      std::swap( ev[0], ev[1] );
      }

    const Precision V = TubenessFunction2D( ev[0], ev[1] );
    vnl_vector<Precision> evn( 2 );

    // adjusting eigenvalues
    // static_cast required to prevent error with gcc 4.1.2
    evn[0]   = 1.0 + ( m_Epsilon - 1.0 ) *
      std::pow( V, static_cast<Precision>( 1.0/m_Sensitivity ) );
    evn[1]   = 1.0 + ( m_Epsilon - 1.0 ) *
      std::pow( V, static_cast<Precision>( 1.0/m_Sensitivity ) );

    vnl_matrix<Precision> LAM( 2, 2 );
    LAM.fill( 0 );
    LAM( 0, 0 ) = evn[0];
    LAM( 1, 1 ) = evn[1];

    const vnl_matrix<Precision> HN = EV * LAM * EV.transpose();

    itxx.Value() = HN( 0, 0 );
    itxy.Value() = HN( 0, 1 );
    ityy.Value() = HN( 1, 1 );
    }
  return;
}


template< class TPixel, unsigned int TDimension >
void
TubeEnhancingDiffusion2DImageFilter<TPixel, TDimension>
::GenerateData( void )
{
  if( m_Verbose )
    {
    std::cout << std::endl <<
      "begin vesselenhancingdiffusion2Dimagefilter ... " << std::endl;
    }

  ProgressReporter progress( this, 0, m_Iterations+4 );

  typedef MinimumMaximumImageFilter<ImageType> MinMaxType;
  typename MinMaxType::Pointer minmax = MinMaxType::New();

  minmax->SetInput( this->GetInput() );

  minmax->Update();

  progress.CompletedPixel();

  const typename ImageType::SpacingType
    ispacing = this->GetInput()->GetSpacing();
  const Precision htmax = 0.5 /
     ( 1.0 / ( ispacing[0] * ispacing[0] )
       + 1.0 / ( ispacing[1] * ispacing[1] ) );

   if( m_TimeStep == NumericTraits<Precision>::Zero )
    {
    m_TimeStep = htmax;
    }

  if( m_TimeStep> htmax )
    {
    std::cerr << "the time step size is too large!" << std::endl;
    this->AllocateOutputs();
    return;
    }

  if( m_Verbose )
    {
    std::cout << "min/max             \t" << minmax->GetMinimum()
              << " " << minmax->GetMaximum() << std::endl;
    std::cout << "iterations/timestep \t" << m_Iterations
              << " " << m_TimeStep << std::endl;
    std::cout << "recalc v            \t" << m_RecalculateTubeness
              << std::endl;
    std::cout << "scales              \t";
    for( unsigned int i=0; i<m_Scales.size(); ++i )
      {
      std::cout << m_Scales[i] << " ";
      }
    std::cout << std::endl;
    std::cout << "eps/omega/sens      \t" << m_Epsilon
              <<  " " << m_Omega << " " << m_Sensitivity << std::endl;
    }

  // cast to precision
  typedef CastImageFilter<ImageType, PrecisionImageType> CT;
  typename CT::Pointer cast = CT::New();
  cast->SetInput( this->GetInput() );
  cast->Update();

  typename PrecisionImageType::Pointer ci = cast->GetOutput();

  progress.CompletedPixel();

  if( m_Verbose )
    {
    std::cout << "start algorithm ... " << std::endl;
    }

  for( m_CurrentIteration=1;
       m_CurrentIteration<=m_Iterations;
       m_CurrentIteration++ )
    {
    VED2DSingleIteration ( ci );
    progress.CompletedPixel();
    }

  typedef MinimumMaximumImageFilter<PrecisionImageType> MMT;
  typename MMT::Pointer mm = MMT::New();
  mm->SetInput( ci );
  mm->Update();

  progress.CompletedPixel();

  if( m_Verbose )
    {
    std::cout << std::endl;
    std::cout << "min/max             \t" << mm->GetMinimum()
              << " " << mm->GetMaximum() << std::endl;
    std::cout << "end vesselenhancingdiffusion2Dimagefilter"
              << std::endl;
    }

  // cast back to pixel type
  this->AllocateOutputs();
  typedef CastImageFilter<PrecisionImageType, ImageType> CTI;
  typename CTI::Pointer casti = CTI::New();
  casti->SetInput( ci );
  casti->GraftOutput( this->GetOutput() );
  casti->Update();
  this->GraftOutput( casti->GetOutput() );

  progress.CompletedPixel();
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeTubeEnhancingDiffusion2DImageFilter_hxx )
