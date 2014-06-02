/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkRegionFromReferenceImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2008-10-13 18:54:28 $
  Version:   $Revision: 1.4 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itktubeRegionFromReferenceImageFilter_txx
#define __itktubeRegionFromReferenceImageFilter_txx

#include "itktubeRegionFromReferenceImageFilter.h"

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
  this->SetExtractionRegion( this->GetReferenceImage()->GetLargestPossibleRegion() );
  Superclass::GenerateOutputInformation();
}
  
/**
 *
 */
template <class TInputImage, class TOutputImage>
void 
RegionFromReferenceImageFilter<TInputImage,TOutputImage>
::SetReferenceImage ( const ReferenceImageType *image )
{
  itkDebugMacro("setting input ReferenceImage to " << image);
  if( image != static_cast<const ReferenceImageType *>(this->GetInput( 1 )) )
    {
    this->ProcessObject::SetNthInput(1, const_cast< ReferenceImageType *>( image ) );
    this->Modified();
    }
}
  
/**
 *
 */
template <class TInputImage, class TOutputImage>
const typename RegionFromReferenceImageFilter<TInputImage,TOutputImage>::ReferenceImageType *
RegionFromReferenceImageFilter<TInputImage,TOutputImage>
::GetReferenceImage() const
{
  Self * surrogate = const_cast< Self * >( this );

  const DataObject * input = surrogate->ProcessObject::GetInput(1);

  const ReferenceImageType * referenceImage = static_cast<const ReferenceImageType *>( input );
  
  return referenceImage;
}

} // end namespace tube

} // end namespace itk

#endif
