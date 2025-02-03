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
#ifndef __tubeComputeTrainingMask_hxx
#define __tubeComputeTrainingMask_hxx


namespace tube
{

/**
 *
 */
template <class TImage, class TLabelMap>
ComputeTrainingMask<TImage, TLabelMap>::ComputeTrainingMask(void)
{
  m_Filter = FilterType::New();
}


/**
 *
 */
template <class TImage, class TLabelMap>
void
ComputeTrainingMask<TImage, TLabelMap>::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  os << indent << "Gap:" << m_Filter->GetGap() << std::endl;
  os << indent << "NotObjectWidth:" << m_Filter->GetNotObjectWidth() << std::endl;
}


} // end namespace tube

#endif
