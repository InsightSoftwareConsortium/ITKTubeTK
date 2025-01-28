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
#ifndef __tubeSegmentTubes_hxx
#define __tubeSegmentTubes_hxx


namespace tube
{
template <class TInputImage>
SegmentTubes<TInputImage>::SegmentTubes(void)
{
  m_Filter = FilterType::New();
  m_RidgeFilter = m_Filter->GetRidgeExtractor();
  m_RadiusFilter = m_Filter->GetRadiusExtractor();

  m_Verbose = false;
  m_Ridgeness = 0;
  m_Intensity = 0;
  m_Roundness = 0;
  m_Curvature = 0;
  m_Levelness = 0;
}

template <class TInputImage>
void
SegmentTubes<TInputImage>::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << m_Filter << std::endl;
}

} // end namespace tube


#endif
