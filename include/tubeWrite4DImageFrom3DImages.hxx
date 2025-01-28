/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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
#ifndef __tubeWrite4DImageFrom3DImages_hxx
#define __tubeWrite4DImageFrom3DImages_hxx


#include "tubeMessage.h"

#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <itksys/SystemTools.hxx>


namespace tube
{

template <class InputImageT>
Write4DImageFrom3DImages<InputImageT>::Write4DImageFrom3DImages()
{
  m_OutputImage = nullptr;

  m_NumberOfInputImages = 0;

  m_FileName = "";
}

template <class InputImageT>
void
Write4DImageFrom3DImages<InputImageT>::SetNumberOfInputImages(unsigned int numInputs)
{
  m_NumberOfInputImages = numInputs;

  m_OutputImage = nullptr;
}

template <class InputImageT>
void
Write4DImageFrom3DImages<InputImageT>::SetNthInputImage(unsigned int outputIndex, const InputImageType * img)
{
  if (m_OutputImage == nullptr)
  {
    m_OutputImage = OutputImageType::New();

    typename InputImageType::RegionType inRegion;
    inRegion = img->GetLargestPossibleRegion();

    typename OutputImageType::SizeType      outSize;
    typename OutputImageType::SpacingType   outSpacing;
    typename OutputImageType::IndexType     outIndex;
    typename OutputImageType::DirectionType outDirection;
    outDirection.SetIdentity();
    for (unsigned int i = 0; i < 3; ++i)
    {
      outSize[i] = inRegion.GetSize()[i];
      outSpacing[i] = img->GetSpacing()[i];
      outIndex[i] = inRegion.GetIndex()[i];
      for (unsigned int j = 0; j < 3; ++j)
      {
        outDirection(i, j) = img->GetDirection()(i, j);
      }
    }
    outSize[3] = m_NumberOfInputImages;
    outSpacing[3] = 1;
    outIndex[3] = 0;
    for (unsigned int j = 0; j < 3; ++j)
    {
      outDirection(3, j) = 0;
    }
    outDirection(3, 3) = 1;

    typename OutputImageType::RegionType outRegion;
    outRegion.SetSize(outSize);
    outRegion.SetIndex(outIndex);
    m_OutputImage->SetRegions(outRegion);
    m_OutputImage->SetSpacing(outSpacing);
    m_OutputImage->SetDirection(outDirection);
    m_OutputImage->Allocate(0);
  }

  itk::ImageRegionConstIterator<InputImageType> inIter(img, img->GetLargestPossibleRegion());
  itk::ImageRegionIterator<OutputImageType>     outIter(m_OutputImage, m_OutputImage->GetLargestPossibleRegion());
  inIter.GoToBegin();
  outIter.GoToBegin();
  unsigned int numPixelsPerSlice = img->GetLargestPossibleRegion().GetSize()[0];
  for (unsigned int i = 1; i < 3; ++i)
  {
    numPixelsPerSlice *= img->GetLargestPossibleRegion().GetSize()[i];
  }
  for (unsigned int i = 0; i < outputIndex * numPixelsPerSlice; ++i)
  {
    ++outIter;
  }
  while (!inIter.IsAtEnd())
  {
    outIter.Set(inIter.Get());
    ++inIter;
    ++outIter;
  }
}

/** Main work happens here */
template <class InputImageT>
void
Write4DImageFrom3DImages<InputImageT>::Update()
{
  typedef itk::ImageFileWriter<OutputImageType> WriterType;
  typename WriterType::Pointer                  writer = WriterType::New();
  writer->SetInput(m_OutputImage);
  writer->SetFileName(m_FileName);
  writer->SetUseCompression(true);
  writer->Update();
}

template <class InputImageT>
void
Write4DImageFrom3DImages<InputImageT>::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Output image = " << m_OutputImage << std::endl;
  os << indent << "File name = " << m_FileName << std::endl;
}

}; // namespace tube

#endif
