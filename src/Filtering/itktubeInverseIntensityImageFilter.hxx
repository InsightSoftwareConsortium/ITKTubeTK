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

#ifndef __itktubeInverseIntensityImageFilter_hxx
#define __itktubeInverseIntensityImageFilter_hxx


#include <itkImageRegionIterator.h>
#include <itkMinimumMaximumImageFilter.h>

namespace itk
{

namespace tube
{


template< class TInputImage >
InverseIntensityImageFilter<TInputImage>
::InverseIntensityImageFilter( void )
{
  m_InverseMaximumIntensity = 0;
}

/** Generate Data */
template< class TInputImage >
void
InverseIntensityImageFilter<TInputImage>
::GenerateData( void )
{
  itkDebugMacro( << "InverseIntensityImageFilter::Generate Data() called." );

  /**************************************/
  /* Get the input and output pointers */
  /************************************/
  InputImageConstPointer  InputImage = this->GetInput();
  RegionType region;             //
  region = InputImage->GetLargestPossibleRegion();    //

  OutputImagePointer OutputImage  = this->GetOutput();
  OutputImage->SetLargestPossibleRegion( region );   //
  OutputImage->SetBufferedRegion( region );         // set the region
  OutputImage->SetRequestedRegion( region );       //
  OutputImage->SetSpacing( InputImage->GetSpacing() );  // set spacing
  OutputImage->SetOrigin( InputImage->GetOrigin() );   //   and origin
  OutputImage->SetDirection( InputImage->GetDirection() );   //   and origin
  OutputImage->Allocate();                          // allocate the image

  InputPixelType myMin;        //
  InputPixelType myMax;       //  declaration variables
  InputPixelType value;      //       calc max/min
  InputPixelType temp;     //

  typename itk::MinimumMaximumImageFilter<InputImageType>::Pointer
    MinMaxFilter;
  MinMaxFilter=itk::MinimumMaximumImageFilter<InputImageType>::New();
  MinMaxFilter->SetInput( InputImage );         //  compute max
  MinMaxFilter->Update();                       // and min
  myMin=MinMaxFilter->GetMinimum();             //  of the input image
  myMax=MinMaxFilter->GetMaximum();             //

  // ** Input max value is given, then use it instead of true max ** //
  if( m_InverseMaximumIntensity ) { myMax = m_InverseMaximumIntensity;  }

  // ** Build Input and output iterators for inversion ** //
  typedef  ImageRegionConstIterator<InputImageType>      InputIteratorType;
  typedef  ImageRegionIterator<OutputImageType>          OutputIteratorType;
  InputIteratorType it_input( InputImage, region );
  OutputIteratorType it_output( OutputImage, region );

  it_input.GoToBegin();
  it_output.GoToBegin();

  // ** Since output is a product of input, only output needs to be tested ** //
  while( !it_output.IsAtEnd() )
    {
    temp = it_input.Get();
    if( temp > myMax )
      {
      value = myMin;
      }
    else
      {
      value = ( myMax + myMin ) - temp;
      }
    OutputPixelType out = ( OutputPixelType ) value;

    it_output.Set( out );
    ++it_output;
    ++it_input;
    }

  itkDebugMacro( << "InverseIntensityImageFilter::Generate Data() finished." );
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeInverseIntensityImageFilter_hxx )
