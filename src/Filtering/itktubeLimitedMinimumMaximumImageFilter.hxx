/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
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
#ifndef itktubeLimitedMinimumMaximumImageFilter_hxx
#define itktubeLimitedMinimumMaximumImageFilter_hxx


#include "itkImageScanlineIterator.h"
#include <mutex>

#include <vector>

namespace itk
{

namespace tube
{

template <typename TInputImage>
LimitedMinimumMaximumImageFilter<TInputImage>::LimitedMinimumMaximumImageFilter()
{
  Self::SetMinimumLimit(NumericTraits<PixelType>::NonpositiveMin());
  Self::SetMaximumLimit(NumericTraits<PixelType>::max());

  Self::SetMinimum(NumericTraits<PixelType>::max());
  Self::SetMaximum(NumericTraits<PixelType>::NonpositiveMin());
}

template <typename TInputImage>
DataObject::Pointer
LimitedMinimumMaximumImageFilter<TInputImage>::MakeOutput(const DataObjectIdentifierType & name)
{
  if (name == "Minimum" || name == "Maximum")
  {
    return PixelObjectType::New();
  }
  return Superclass::MakeOutput(name);
}

template <typename TInputImage>
void
LimitedMinimumMaximumImageFilter<TInputImage>::BeforeStreamedGenerateData()
{
  Superclass::BeforeStreamedGenerateData();

  m_ThreadMin = m_MaximumLimit;
  m_ThreadMax = m_MinimumLimit;
}

template <typename TInputImage>
void
LimitedMinimumMaximumImageFilter<TInputImage>::AfterStreamedGenerateData()
{
  Superclass::AfterStreamedGenerateData();

  this->SetMinimum(m_ThreadMin);
  this->SetMaximum(m_ThreadMax);
}

template <typename TInputImage>
void
LimitedMinimumMaximumImageFilter<TInputImage>::ThreadedStreamedGenerateData(const RegionType & regionForThread)
{
  if (regionForThread.GetNumberOfPixels() == 0)
  {
    return;
  }

  PixelType localMin = m_MaximumLimit;
  PixelType localMax = m_MinimumLimit;

  ImageScanlineConstIterator<TInputImage> it(this->GetInput(), regionForThread);


  // do the work
  while (!it.IsAtEnd())
  {
    // Handle the odd pixel separately
    if (regionForThread.GetSize(0) % 2 == 1)
    {
      const PixelType value = it.Get();
      localMin = std::min(value, localMin);
      localMax = std::max(value, localMax);
      ++it;
    }

    while (!it.IsAtEndOfLine())
    {
      const PixelType value1 = it.Get();
      ++it;
      const PixelType value2 = it.Get();
      ++it;

      if (value1 > value2)
      {
        localMax = std::max(value1, localMax);
        localMin = std::min(value2, localMin);
      }
      else
      {
        localMax = std::max(value2, localMax);
        localMin = std::min(value1, localMin);
      }
    }
    it.NextLine();
  }

  std::lock_guard<std::mutex> mutexHolder(m_Mutex);
  m_ThreadMin = std::min(localMin, m_ThreadMin);
  m_ThreadMax = std::max(localMax, m_ThreadMax);
}

template <typename TImage>
void
LimitedMinimumMaximumImageFilter<TImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Minimum Limit: " << static_cast<typename NumericTraits<PixelType>::PrintType>(m_MinimumLimit)
     << std::endl;
  os << indent << "Maximum Limit: " << static_cast<typename NumericTraits<PixelType>::PrintType>(m_MaximumLimit)
     << std::endl;
  os << indent << "Minimum: " << static_cast<typename NumericTraits<PixelType>::PrintType>(this->GetMinimum())
     << std::endl;
  os << indent << "Maximum: " << static_cast<typename NumericTraits<PixelType>::PrintType>(this->GetMaximum())
     << std::endl;
}
} // end namespace tube
} // end namespace itk
#endif
