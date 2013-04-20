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
#ifndef __itkTubeSupervisedClassifierFilter_txx
#define __itkTubeSupervisedClassifierFilter_txx

#include <limits>
#include <iostream>

#include "itkTubeSupervisedClassifierFilter.h"

#include "itkTimeProbesCollectorBase.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "tubeMatrixMath.h"

namespace itk
{

namespace tube
{

template< class ImageT, class LabelmapT >
SupervisedClassifierFilter< ImageT, LabelmapT >
::SupervisedClassifierFilter()
{
  m_Labelmap = NULL;

  m_FeatureVectorGenerator = NULL;

  m_ObjectIdList.clear();

  m_ProgressProcessInfo = NULL;
  m_ProgressFraction = 1.0;
  m_ProgressStart = 0;
}

template< class ImageT, class LabelmapT >
SupervisedClassifierFilter< ImageT, LabelmapT >
::~SupervisedClassifierFilter()
{
}

template < class ImageT, class LabelmapT >
void
SupervisedClassifierFilter< ImageT, LabelmapT >
::SetFeatureVectorGenerator( FeatureVectorGenerator::Pointer fGen )
{
  m_FeatureVectorGenerator = fGen;
}

template < class ImageT, class LabelmapT >
typename SupervisedClassifierFilter< ImageT, LabelmapT >::FeatureVectorGeneratorType::Pointer *
SupervisedClassifierFilter< ImageT, LabelmapT >
::GetFeatureVectorGenerator( void )
{
  return m_FeatureVectorGenerator;
}

template < class ImageT, class LabelmapT >
void
SupervisedClassifierFilter< ImageT, LabelmapT >
::SetObjectId( ObjectIdType objectId )
{
  m_ObjectIdList.clear();
  m_ObjectIdList.push_back( objectId );
}

template < class ImageT, class LabelmapT >
void
SupervisedClassifierFilter< ImageT, LabelmapT >
::AddObjectId( ObjectIdType objectId )
{
  m_ObjectIdList.push_back( objectId );
}

template < class ImageT, class LabelmapT >
unsigned int
SupervisedClassifierFilter< ImageT, LabelmapT >
::GetNumberOfObjectIds( void )
{
   return m_ObjectIdList.size();
}

template < class ImageT, class LabelmapT >
typename SupervisedClassifierFilter< ImageT, LabelmapT >::ObjectIdType
SupervisedClassifierFilter< ImageT, LabelmapT >
::GetObjectId( unsigned int num )
{
  if( num < m_ObjectIdList.size() )
    {
    return m_ObjectIdList[ num ];
    }
  else
    {
    // voidId
    return m_ObjectIdList[ m_ObjectIdList.size()-1 ];
    }
}

template < class ImageT, class LabelmapT >
void
SupervisedClassifierFilter< ImageT, LabelmapT >
::SetProgressProcessInformation( void * processInfo, double fraction,
  double start )
{
  m_ProgressProcessInfo = processInfo;
  m_ProgressFraction = fraction;
  m_ProgressStart = start;
}


template <class ImageT, class LabelmapT >
void
SupervisedClassifierFilter< ImageT, LabelmapT >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_Labelmap )
    {
    os << indent << "Labelmap = " << m_Labelmap << std::endl;
    }
  else
    {
    os << indent << "Labelmap = NULL" << std::endl;
    }

  os << indent << "ObjectIdList.size = " << m_ObjectIdList.size()
    << std::endl;

  if( m_ProgressProcessInfo != NULL )
    {
    os << indent << "ProgressProcessInfo = Set" << std::endl;
    }
  else
    {
    os << indent << "ProgressProcessInfo = NULL" << std::endl;
    }
  os << indent << "ProgressFraction = " << m_ProgressFraction << std::endl;
  os << indent << "ProgressStart = " << m_ProgressStart << std::endl;

}

//----------------------------------------------------------------------------
template <class ImageT, class LabelmapT >
void
SupervisedClassifierFilter< ImageT, LabelmapT >
::BeforeThreadedGenerateData()
{
  /*
  Example multi-threaded summation code

  int numberOfThreads = this->GetNumberOfThreads();

  // Initialize statistics about difference image.
  m_MeanDifference = NumericTraits<RealType>::Zero;
  m_TotalDifference = NumericTraits<AccumulateType>::Zero;
  m_NumberOfPixelsWithDifferences = 0;

  // Resize the thread temporaries
  m_ThreadDifferenceSum.SetSize(numberOfThreads);
  m_ThreadNumberOfPixelsWithDifferences.SetSize(numberOfThreads);

  // Initialize the temporaries
  m_ThreadDifferenceSum.Fill(NumericTraits<AccumulateType>::Zero);
  m_ThreadNumberOfPixelsWithDifferences.Fill(0);
  */
}

//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage>
void
DifferenceImageFilter2<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType &threadRegion,
  ThreadIdType threadId)
{
  /*
  Example multi-threaded summation code
  typedef ConstNeighborhoodIterator<InputImageType>   SmartIterator;
  typedef ImageRegionConstIterator<InputImageType>    InputIterator;
  typedef ImageRegionIterator<OutputImageType>        OutputIterator;
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>
                                                      FacesCalculator;
  typedef typename FacesCalculator::RadiusType        RadiusType;
  typedef typename FacesCalculator::FaceListType      FaceListType;
  typedef typename FaceListType::iterator             FaceListIterator;
  typedef typename InputImageType::PixelType          InputPixelType;

  // Prepare standard boundary condition.
  ZeroFluxNeumannBoundaryCondition<InputImageType> nbc;

  // Get a pointer to each image.
  const InputImageType* validImage = this->GetInput(0);
  const InputImageType* testImage = this->GetInput(1);
  OutputImageType* outputPtr = this->GetOutput();

  // Create a radius of pixels.
  RadiusType radius;
  if(m_ToleranceRadius > 0)
    {
    radius.Fill(m_ToleranceRadius);
    }
  else
    {
    radius.Fill(0);
    }

  // Find the data-set boundary faces.
  FacesCalculator boundaryCalculator;
  FaceListType faceList = boundaryCalculator(testImage, threadRegion, radius);

  // Support progress methods/callbacks.
  ProgressReporter progress(this, threadId, threadRegion.GetNumberOfPixels());

  // Process the internal face and each of the boundary faces.
  for (FaceListIterator face = faceList.begin(); face != faceList.end(); ++face)
    {
    SmartIterator test(radius, testImage, *face); // Iterate over test image.
    InputIterator valid(validImage, *face);       // Iterate over valid image.
    OutputIterator out(outputPtr, *face);         // Iterate over output image.
    if( !test.GetNeedToUseBoundaryCondition() || !m_IgnoreBoundaryPixels )
      {
      test.OverrideBoundaryCondition(&nbc);

      for(valid.GoToBegin(), test.GoToBegin(), out.GoToBegin();
          !valid.IsAtEnd();
          ++valid, ++test, ++out)
        {
        // Get the current valid pixel.
        InputPixelType t = valid.Get();

        //  Assume a good match - so test center pixel first, for speed
        RealType difference = static_cast<RealType>(t) - test.GetCenterPixel();
        RealType absDifference = difference;
        if(NumericTraits<RealType>::IsNegative(difference))
          {
          absDifference = -absDifference;
          }
        RealType minimumDifference = difference;
        RealType minimumAbsDifference = absDifference;

        unsigned int neighborhoodSize = test.Size();
        for (unsigned int i=0; i < neighborhoodSize; ++i)
          {
          // Use the RealType for the difference to make sure we get the
          // sign.
          difference = static_cast<RealType>(t) - test.GetPixel(i);

          absDifference = difference;
          if(NumericTraits<RealType>::IsNegative(difference))
            {
            absDifference = -absDifference;
            }
          if(absDifference < minimumAbsDifference)
            {
            minimumDifference = difference;
            minimumAbsDifference = absDifference;
            }
          }

        // Store the minimum difference value in the output image.
        out.Set( static_cast<OutputPixelType>(minimumDifference) );

        m_ThreadDifferenceSum[threadId] += minimumAbsDifference;
        if(minimumAbsDifference > m_DifferenceThreshold)
          {
          m_ThreadNumberOfPixelsWithDifferences[threadId]++;
          }

        // Update progress.
        progress.CompletedPixel();
        }
      }
    else
      {
      for(out.GoToBegin(); !out.IsAtEnd(); ++out)
        {
        out.Set(NumericTraits<OutputPixelType>::Zero);
        progress.CompletedPixel();
        }
      }
    }
  */
}

//----------------------------------------------------------------------------
template <class TInputImage, class TOutputImage>
void
DifferenceImageFilter2<TInputImage, TOutputImage>
::AfterThreadedGenerateData()
{
  /*
  Example multi-threaded summation code

  // Set statistics about difference image.
  int numberOfThreads = this->GetNumberOfThreads();
  for(int i=0; i < numberOfThreads; ++i)
    {
    m_TotalDifference += m_ThreadDifferenceSum[i];
    m_NumberOfPixelsWithDifferences +=
      m_ThreadNumberOfPixelsWithDifferences[i];
    }

  // Get the total number of pixels processed in the region.
  // This is different from the m_TotalNumberOfPixels which
  // is the number of pixels that actually have differences
  // above the intensity threshold.
  OutputImageRegionType region = this->GetOutput()->GetRequestedRegion();
  AccumulateType numberOfPixels = region.GetNumberOfPixels();

  // Calculate the mean difference.
  m_MeanDifference = m_TotalDifference / numberOfPixels;
  */
}

}

}

#endif //SupervisedClassifierFilter_txx
