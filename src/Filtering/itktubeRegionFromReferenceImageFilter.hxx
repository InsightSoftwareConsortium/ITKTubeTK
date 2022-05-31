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
#ifndef __itktubeRegionFromReferenceImageFilter_hxx
#define __itktubeRegionFromReferenceImageFilter_hxx


namespace itk {

namespace tube {

template <class TInputImage, class TOutputImage>
void
RegionFromReferenceImageFilter<TInputImage, TOutputImage>
::GenerateOutputInformation()
{
  if( !this->GetInput() || !this->GetReferenceImage() )
    {
    return;
    }

  // Superclass::Superclass::GenerateOutputInformation();
  this->SetExtractionRegion(
    this->GetReferenceImage()->GetLargestPossibleRegion() );
  Superclass::GenerateOutputInformation();
}
 
/**
 *
 */
template <class TInputImage, class TOutputImage>
void
RegionFromReferenceImageFilter<TInputImage, TOutputImage>
::SetReferenceImage ( const ReferenceImageType *image )
{
  itkDebugMacro( "setting input ReferenceImage to " << image );
  if( image != static_cast<const ReferenceImageType *>(
    this->GetInput( 1 ) ) )
    {
    this->ProcessObject::SetNthInput( 1,
      const_cast< ReferenceImageType *>( image ) );
    this->Modified();
    }
}
 
/**
 *
 */
template <class TInputImage, class TOutputImage>
const typename
RegionFromReferenceImageFilter<TInputImage, TOutputImage>::
ReferenceImageType *
RegionFromReferenceImageFilter<TInputImage, TOutputImage>
::GetReferenceImage() const
{
  Self * surrogate = const_cast< Self * >( this );

  const DataObject * input = surrogate->ProcessObject::GetInput( 1 );

  const ReferenceImageType * referenceImage =
    static_cast<const ReferenceImageType *>( input );
 
  return referenceImage;
}

} // end namespace tube

} // end namespace itk

#endif
