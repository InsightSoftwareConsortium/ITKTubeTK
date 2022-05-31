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
#ifndef __itktubeGaussianDerivativeFilter_hxx
#define __itktubeGaussianDerivativeFilter_hxx

namespace itk
{
namespace tube
{
//----------------------------------------------------------------------------

template< typename TInputImage, typename TOutputImage >
GaussianDerivativeFilter<TInputImage, TOutputImage>
::GaussianDerivativeFilter()
{
  this->m_Orders.Fill ( 0 );
  this->m_Sigmas.Fill ( 0 );
}

template< typename TInputImage, typename TOutputImage >
void
GaussianDerivativeFilter<TInputImage, TOutputImage>
::PrintSelf ( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf ( os, indent );

  os << indent << "Orders              : " << m_Orders << std::endl;
  os << indent << "Sigmas               : " << m_Sigmas << std::endl;
}

} // End namespace tube
} // End namespace itk

#endif
