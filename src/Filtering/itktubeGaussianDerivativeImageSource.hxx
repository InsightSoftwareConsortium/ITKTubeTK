/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/
#ifndef __itktubeGaussianDerivativeImageSource_hxx
#define __itktubeGaussianDerivativeImageSource_hxx

#include <itkGaussianSpatialFunction.h>
#include <itkImageRegionIterator.h>
#include <itkObjectFactory.h>
#include <itkProgressReporter.h>

namespace itk
{

namespace tube
{

//----------------------------------------------------------------------------
template< typename TOutputImage >
GaussianDerivativeImageSource< TOutputImage >
::GaussianDerivativeImageSource()
{
  // Gaussian parameters, defined so that the gaussian
  // is centered in the default image
  m_Index.Fill( 0 );

  m_Mean.Fill( 0 );
  m_Sigmas.Fill( 1 );
  m_Orders.Fill( 0 );
}

//----------------------------------------------------------------------------
template< typename TOutputImage >
void
GaussianDerivativeImageSource< TOutputImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Index: " << m_Index << std::endl;

  os << indent << "Gaussian order: " << m_Orders << std::endl;

  os << indent << "Gaussian sigma: " << m_Sigmas << std::endl;

  os << indent << "Gaussian mean: " << m_Mean << std::endl;
}

//----------------------------------------------------------------------------
template< typename TOutputImage >
void
GaussianDerivativeImageSource< TOutputImage >
::SetParameters( const ParametersType & parameters )
{
  OrdersType orders;
  SigmasType sigmas;
  PointType mean;
  for( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
    {
    orders[i] = parameters[i];
    sigmas[i]  = parameters[i + TOutputImage::ImageDimension];
    mean[i]  = parameters[i + 2*TOutputImage::ImageDimension];
    }
  this->SetOrders( orders );
  this->SetSigmas( sigmas );
  this->SetMean( mean );
}

//----------------------------------------------------------------------------
template< typename TOutputImage >
typename GaussianDerivativeImageSource< TOutputImage >::ParametersType
GaussianDerivativeImageSource< TOutputImage >
::GetParameters() const
{
  ParametersType parameters( 2*TOutputImage::ImageDimension );
  for( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
    {
    parameters[i] = m_Orders[i];
    parameters[i + TOutputImage::ImageDimension] = m_Sigmas[i];
    parameters[i + 2*TOutputImage::ImageDimension] = m_Mean[i];
    }

  return parameters;
}

//----------------------------------------------------------------------------
template< typename TOutputImage >
unsigned int
GaussianDerivativeImageSource< TOutputImage >
::GetNumberOfParameters() const
{
  return 3*TOutputImage::ImageDimension;
}

//----------------------------------------------------------------------------
template< typename TOutputImage >
void
GaussianDerivativeImageSource< TOutputImage >
::GenerateOutputInformation()
{
  OutputImageType *output = this->GetOutput( 0 );

  typename OutputImageType::RegionType largestPossibleRegion;
  largestPossibleRegion.SetSize( this->GetSize() );
  largestPossibleRegion.SetIndex( this->m_Index );
  output->SetLargestPossibleRegion( largestPossibleRegion );

  output->SetSpacing( this->GetSpacing() );
  output->SetOrigin( this->GetOrigin() );
  output->SetDirection( this->GetDirection() );
}

//----------------------------------------------------------------------------
template< typename TOutputImage >
void
GaussianDerivativeImageSource< TOutputImage >
::GenerateData()
{
  TOutputImage * outputPtr = this->GetOutput();
  // allocate the output buffer
  outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
  outputPtr->Allocate();

  // Create an iterator that will walk the output region
  typedef itk::ImageRegionIterator< TOutputImage > OutputIterator;
  OutputIterator outIt = OutputIterator( outputPtr,
                                         outputPtr->GetRequestedRegion() );

  // The position at which the function is evaluated
  Point< double, TOutputImage::ImageDimension > evalPoint;

  ProgressReporter progress( this, 0,
                             outputPtr->GetRequestedRegion()
                             .GetNumberOfPixels() );
  double prefixDenom = 1.0;
  const double squareRootOfTwoPi = std::sqrt( 2.0 * vnl_math::pi );
  for( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
    {
    prefixDenom *= m_Sigmas[i] * squareRootOfTwoPi;
    }
  double initPrefixDenom = prefixDenom;

  double total = 0;
  while( !outIt.IsAtEnd() )
    {
    typename TOutputImage::IndexType index = outIt.GetIndex();
    outputPtr->TransformIndexToPhysicalPoint( index, evalPoint );

    prefixDenom = initPrefixDenom;

    for( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
      {
      if( m_Orders[i] != 0 )
        {
        prefixDenom *= std::pow( m_Sigmas[i], 2*m_Orders[i] )
          / ( std::pow( ( -( evalPoint[i] - m_Mean[i] ) ), m_Orders[i] ) -
          ( m_Orders[i] == 2 ? std::pow( m_Sigmas[1], m_Orders[i] ) : 0 ) );
        }
      }
    double suffixExp = 0;
    for( unsigned int i = 0; i < TOutputImage::ImageDimension; i++ )
      {
      suffixExp += ( evalPoint[i] - m_Mean[i] )
                   * ( evalPoint[i] - m_Mean[i] )
                   / ( 2 * m_Sigmas[i] * m_Sigmas[i] );
      }

    double value = ( 1 / prefixDenom ) * std::exp( -suffixExp );
    total += std::abs( value );

    // Set the pixel value to the function value
    outIt.Set( ( typename TOutputImage::PixelType )value );

    progress.CompletedPixel();

    ++outIt;
    }

  outIt.GoToBegin();
  while( !outIt.IsAtEnd() )
    {
    outIt.Set( outIt.Get() / total );
    ++outIt;
    }
}
} // End namespace tube

} // End namespace itk

#endif
