/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/

#ifndef __itktubeArrayFireGlueUtilities_h
#define __itktubeArrayFireGlueUtilities_h

#include "itkImportImageFilter.h"
#include "itkExceptionObject.h"

#include <arrayfire.h>

namespace itk
{

template <typename ImageType>
void convertITKImageToArrayFire ( const ImageType * pItkImage,
  af::array &afArr )
{
  typename ImageType::SizeType imageSize =
    pItkImage->GetLargestPossibleRegion().GetSize();

  if( ImageType::ImageDimension == 1 )
    {
    afArr = af::array ( imageSize[0], pItkImage->GetBufferPointer() );
    }
  else if( ImageType::ImageDimension == 2 )
    {
    afArr = af::array ( imageSize[0], imageSize[1],
      pItkImage->GetBufferPointer() );
    }
  else if( ImageType::ImageDimension == 3 )
    {
    afArr = af::array ( imageSize[0], imageSize[1], imageSize[2],
      pItkImage->GetBufferPointer() );
    }
  else
    {
    itk::InvalidArgumentError e ( __FILE__, __LINE__ );
    e.SetDescription ( "Only Dimensions upto 3 is not supported" );
    e.SetLocation ( "convertITKImageToArrayFire" );
    throw e;
    }
}

template <typename ImageType>
void convertArrayFireImageToITK ( af::array &afArr,
  typename ImageType::Pointer &pItkImage,
  const ImageType * pItkReferenceImage )
{
  typedef typename ImageType::PixelType PixelType;
  const unsigned int Dimension = ImageType::ImageDimension;

  typedef itk::ImportImageFilter< PixelType, Dimension > ImportFilterType;
  typename ImportFilterType::Pointer pImportFilter = ImportFilterType::New();

  pImportFilter->SetRegion ( pItkReferenceImage->GetLargestPossibleRegion() );
  pImportFilter->SetOrigin ( pItkReferenceImage->GetOrigin() );
  pImportFilter->SetSpacing ( pItkReferenceImage->GetSpacing() );
  pImportFilter->SetDirection ( pItkReferenceImage->GetDirection() );

  PixelType *pArrayFireImageBuffer = afArr.host<PixelType>();
  pImportFilter->SetImportPointer ( pArrayFireImageBuffer, afArr.elements(),
                                    false );
  pImportFilter->Update();

  pItkImage = pImportFilter->GetOutput();
}

} // End namespace itk

#endif
