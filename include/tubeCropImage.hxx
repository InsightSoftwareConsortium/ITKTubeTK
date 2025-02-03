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
#ifndef __tubeCropImage_hxx
#define __tubeCropImage_hxx


namespace tube
{
template <typename TInputImage, typename TOutputImage>
CropImage<TInputImage, TOutputImage>::CropImage(void)
{
  m_Filter = FilterType::New();
}

template <typename TInputImage, typename TOutputImage>
void
CropImage<TInputImage, TOutputImage>::SetSplitInput(InputIndexType splitIndex, InputIndexType roiIndex)
{
  m_Filter->SetSplitInput(splitIndex, roiIndex);
}

template <typename TInputImage, typename TOutputImage>
void
CropImage<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << m_Filter << std::endl;
}

} // end namespace tube


#endif
