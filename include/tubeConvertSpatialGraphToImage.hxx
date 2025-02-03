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
#ifndef __tubeConvertSpatialGraphToImage_hxx
#define __tubeConvertSpatialGraphToImage_hxx


namespace tube
{
template <typename TInputImage, typename TOutputImage>
ConvertSpatialGraphToImage<TInputImage, TOutputImage>::ConvertSpatialGraphToImage(void)
{
  m_Filter = FilterType::New();
}

template <typename TInputImage, typename TOutputImage>
void
ConvertSpatialGraphToImage<TInputImage, TOutputImage>::SetAdjacencyMatrix(vnl_matrix<double> a)
{
  m_Filter->SetAdjacencyMatrix(a);
}

template <typename TInputImage, typename TOutputImage>
void
ConvertSpatialGraphToImage<TInputImage, TOutputImage>::SetBranchnessVector(vnl_vector<double> b)
{
  m_Filter->SetBranchnessVector(b);
}

template <typename TInputImage, typename TOutputImage>
void
ConvertSpatialGraphToImage<TInputImage, TOutputImage>::SetRadiusVector(vnl_vector<double> r)
{
  m_Filter->SetRadiusVector(r);
}

template <typename TInputImage, typename TOutputImage>
void
ConvertSpatialGraphToImage<TInputImage, TOutputImage>::SetCentralityVector(vnl_vector<double> c)
{
  m_Filter->SetCentralityVector(c);
}

template <typename TInputImage, typename TOutputImage>
void
ConvertSpatialGraphToImage<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << m_Filter << std::endl;
}

} // end namespace tube

#endif
