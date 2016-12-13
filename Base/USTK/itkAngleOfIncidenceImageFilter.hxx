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

#ifndef __itkAngleOfIncidenceImageFilter_hxx
#define __itkAngleOfIncidenceImageFilter_hxx

#include "itkAngleOfIncidenceImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkProgressReporter.h>
#include <itkResampleImageFilter.h>

namespace itk
{
/**
 * Constructor
 */
template< class TInputImage, class TOutputImage >
AngleOfIncidenceImageFilter< TInputImage, TOutputImage >
::AngleOfIncidenceImageFilter( void )
{
  m_UltrasoundProbeOrigin.Fill( 0 );

  //eigenvector analysis filter
  m_EigenVectorAnalysisFilter = EigenVectorAnalysisFilterType::New();
  m_EigenVectorAnalysisFilter->SetDimension( ImageDimension );
  m_EigenVectorAnalysisFilter->OrderEigenValuesBy(
      EigenVectorAnalysisFilterType::FunctorType::OrderByValue );

  //eigenvalue analysis filter
  m_EigenValueAnalysisFilter = EigenValueAnalysisFilterType::New();
  m_EigenValueAnalysisFilter->SetDimension( ImageDimension );
  m_EigenValueAnalysisFilter->OrderEigenValuesBy(
        EigenValueAnalysisFilterType::FunctorType::OrderByValue );

  //Hessian filter
  m_HessianFilter = HessianFilterType::New();

  //eigenvector image
  m_PrimaryEigenVectorImage = EigenVectorImageType::New();
  unsigned int vectorLength = 3; // Eigenvector length
  m_PrimaryEigenVectorImage->SetVectorLength ( vectorLength );

}

template< class TInputImage, class TOutputImage >
void
AngleOfIncidenceImageFilter< TInputImage, TOutputImage >
::GenerateData( void )
{
  typedef typename TInputImage::PointType         InputImagePointType;
  typedef typename TInputImage::SpacingType       InputImageSpacingType;
  typedef ImageRegionConstIterator< TInputImage > InputImageConstIteratorType;

  typedef ImageRegionIterator< TOutputImage > OutputImageIteratorType;


  typename TInputImage::ConstPointer inputImagePtr(
    dynamic_cast< const TInputImage  * >(
      this->ProcessObject::GetInput( 0 ) ) );
  typename TOutputImage::Pointer outputImagePtr(
    dynamic_cast< TOutputImage * >(
      this->ProcessObject::GetOutput( 0 ) ) );

  Point< double> originInput;
  originInput.Fill( 0.0 );
  outputImagePtr->SetOrigin( originInput );
  outputImagePtr->SetSpacing( inputImagePtr->GetSpacing() );
  outputImagePtr->SetDirection( inputImagePtr->GetDirection() );

  outputImagePtr->SetRequestedRegion( inputImagePtr->GetRequestedRegion() );
  outputImagePtr->SetBufferedRegion( inputImagePtr->GetBufferedRegion() );
  outputImagePtr->SetLargestPossibleRegion( inputImagePtr->GetLargestPossibleRegion() );
  outputImagePtr->Allocate();
  outputImagePtr->FillBuffer( 0 );

  InputImageConstIteratorType inputIt( inputImagePtr,
    inputImagePtr->GetLargestPossibleRegion() );
  OutputImageIteratorType     outputIt( outputImagePtr,
    outputImagePtr->GetLargestPossibleRegion() );

  outputIt.GoToBegin();
  inputIt.GoToBegin();

  //Compute the Normal Vector Image
  this->ComputeNormalVectorImage();

  //Iterator for the primary eigenvector image with the largest eigenvector
  itk::ImageRegionIterator<EigenVectorImageType> primaryEigenVectorImageIterator;
  primaryEigenVectorImageIterator = itk::ImageRegionIterator<EigenVectorImageType>(
      m_PrimaryEigenVectorImage, m_PrimaryEigenVectorImage->GetRequestedRegion() );
  primaryEigenVectorImageIterator.GoToBegin();


  InputImageSpacingType inputImageSpacing;
  InputImagePointType   inputImageOrigin;

  while( !inputIt.IsAtEnd() )
    {
    //Compute the angle and set it to output iterator.

    itk::VariableLengthVector<double> vectorPixel;
    vectorPixel.SetSize( 3 );
    vectorPixel[0] = primaryEigenVectorImageIterator.Get()[0];
    vectorPixel[1] = primaryEigenVectorImageIterator.Get()[1];
    vectorPixel[2] = primaryEigenVectorImageIterator.Get()[2];

    itk::Vector<double,3>  primaryEigenVector;

    primaryEigenVector[0] = vectorPixel[0];
    primaryEigenVector[1] = vectorPixel[1];
    primaryEigenVector[2] = vectorPixel[2];

    //Vector from the probe origin to the surface voxel
    itk::Vector<double, 3>    beamVector;

    beamVector[0] = inputIt.GetIndex()[0] - m_UltrasoundProbeOrigin[0];
    beamVector[1] = inputIt.GetIndex()[1] - m_UltrasoundProbeOrigin[1];
    beamVector[2] = inputIt.GetIndex()[2] - m_UltrasoundProbeOrigin[2];

    /*
    std::cout << "Beam vector for( "  << inputIt.GetIndex()[0]  << ","
                                     << inputIt.GetIndex()[1]  << ","
                                     << inputIt.GetIndex()[2]  << " ):=( "
                                     << beamVector[0]        << ","
                                     << beamVector[1]        << ","
                                     << beamVector[2]        << " )" << std::endl;
    */

    //Normalize the vectors
    const double zeroNormVectorTolerance = 1e-8;

    if( ( beamVector.GetNorm() > zeroNormVectorTolerance ) &&
         ( primaryEigenVector.GetNorm() > zeroNormVectorTolerance ) )
      {
      beamVector.Normalize();
      primaryEigenVector.Normalize();

      //compute dot product
      double dotProduct = beamVector*primaryEigenVector;

      outputIt.Set( vnl_math_abs( dotProduct ) );
      }

    ++inputIt;
    ++outputIt;
    ++primaryEigenVectorImageIterator;
    }
}

template< class TInputImage, class TOutputImage >
void
AngleOfIncidenceImageFilter< TInputImage, TOutputImage >
::ComputeNormalVectorImage( void )
{
  m_HessianFilter->SetInput( this->GetInput() );

  double sigma = 0.5;
  m_HessianFilter->SetSigma ( sigma );
  m_HessianFilter->Update();

  m_EigenValueAnalysisFilter->SetInput( m_HessianFilter->GetOutput() );
  m_EigenValueAnalysisFilter->Update();

  m_EigenVectorAnalysisFilter->SetInput( m_HessianFilter->GetOutput() );
  m_EigenVectorAnalysisFilter->Update();

  //Generate an image with eigenvector pixel that correspond to the largest eigenvalue
  typename EigenVectorMatrixImageType::ConstPointer eigenVectorImage =
                    m_EigenVectorAnalysisFilter->GetOutput();

  typename EigenVectorMatrixImageType::RegionType region;
  region.SetSize( eigenVectorImage->GetLargestPossibleRegion().GetSize() );
  region.SetIndex( eigenVectorImage->GetLargestPossibleRegion().GetIndex() );
  m_PrimaryEigenVectorImage->SetRegions( region );
  m_PrimaryEigenVectorImage->SetOrigin( eigenVectorImage->GetOrigin() );
  m_PrimaryEigenVectorImage->SetSpacing( eigenVectorImage->GetSpacing() );
  m_PrimaryEigenVectorImage->Allocate();

  //Fill up the buffer with null vector
  unsigned int vectorLength = 3;
  itk::VariableLengthVector< double > nullVector( vectorLength );
  for( unsigned int i=0; i < vectorLength; i++ )
    {
    nullVector[i] = 0.0;
    }
  m_PrimaryEigenVectorImage->FillBuffer( nullVector );

  //Setup the iterators
  //
  //Iterator for the eigenvector matrix image
  itk::ImageRegionConstIterator<EigenVectorMatrixImageType> eigenVectorImageIterator;
  eigenVectorImageIterator = itk::ImageRegionConstIterator<EigenVectorMatrixImageType>(
      eigenVectorImage, eigenVectorImage->GetRequestedRegion() );
  eigenVectorImageIterator.GoToBegin();

  //Iterator for the output image with the largest eigenvector
  itk::ImageRegionIterator<EigenVectorImageType> primaryEigenVectorImageIterator;
  primaryEigenVectorImageIterator = itk::ImageRegionIterator<EigenVectorImageType>(
      m_PrimaryEigenVectorImage, m_PrimaryEigenVectorImage->GetRequestedRegion() );
  primaryEigenVectorImageIterator.GoToBegin();

  //Iterator for the eigenvalue image
  typename EigenValueImageType::ConstPointer eigenImage
    = m_EigenValueAnalysisFilter->GetOutput();
  itk::ImageRegionConstIterator<EigenValueImageType> eigenValueImageIterator;
  eigenValueImageIterator = itk::ImageRegionConstIterator<EigenValueImageType>(
      eigenImage, eigenImage->GetRequestedRegion() );
  eigenValueImageIterator.GoToBegin();

  const double toleranceEigenValues = 1e-4;

  while( !eigenValueImageIterator.IsAtEnd() )
    {
    // Get the eigenvalue
    EigenValueArrayType eigenValue;
    eigenValue = eigenValueImageIterator.Get();

    // Find the largest eigenvalue
    double largest = vnl_math_abs( eigenValue[0] );
    unsigned int largestEigenValueIndex=0;

    for( unsigned int i=1; i <=2; i++ )
      {
      if(  vnl_math_abs( eigenValue[i] > largest ) )
        {
        largest = vnl_math_abs( eigenValue[i] );
        largestEigenValueIndex = i;
        }
      }

    EigenVectorMatrixType   matrixPixel;
    matrixPixel = eigenVectorImageIterator.Get();

    /*
    std::cout << "EigenValues( " << eigenValueImageIterator.GetIndex()[0] << ","
                                         << eigenValueImageIterator.GetIndex()[1] << ","
                                         << eigenValueImageIterator.GetIndex()[2] <<" )\t="
                                         << smallest << ","
                                         << largest << " )" << std::endl;
    */
    if( vnl_math_abs( largest ) >  toleranceEigenValues  )
      {
      //Assuming eigenvectors are rows
      itk::VariableLengthVector<double> primaryEigenVector( vectorLength );
      for( unsigned int i=0; i < vectorLength; i++ )
        {
        primaryEigenVector[i] = matrixPixel[largestEigenValueIndex][i];
        }

      primaryEigenVectorImageIterator.Set( primaryEigenVector );
      }

    ++eigenValueImageIterator;
    ++eigenVectorImageIterator;
    ++primaryEigenVectorImageIterator;
    }
}


template< class TInputImage, class TOutputImage >
void
AngleOfIncidenceImageFilter< TInputImage, TOutputImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Ultrasound origin vector : "
     << static_cast< typename NumericTraits< VectorType >::PrintType >
  ( m_UltrasoundProbeOrigin )
  << std::endl;
}

} // End namespace itk

#endif // End !defined( __itkAngleOfIncidenceImageFilter_hxx )
