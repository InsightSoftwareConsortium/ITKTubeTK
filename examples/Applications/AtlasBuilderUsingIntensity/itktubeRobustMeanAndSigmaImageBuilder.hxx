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

#ifndef __itktubeRobustMeanAndSigmaImageBuilder_hxx
#define __itktubeRobustMeanAndSigmaImageBuilder_hxx


template< class TInputImageType, class TOutputMeanImageType,
  class TOutputSigmaImageType >
RobustMeanAndSigmaImageBuilder< TInputImageType, TOutputMeanImageType,
  TOutputSigmaImageType >
::RobustMeanAndSigmaImageBuilder( void )
: m_NumberOfOutlierImagesToRemove( 0 ),
  m_TotalNumberOfImages( 0 )
{
}

template< class TInputImageType, class TOutputMeanImageType,
  class TOutputSigmaImageType >
void
RobustMeanAndSigmaImageBuilder< TInputImageType, TOutputMeanImageType,
  TOutputSigmaImageType >
::BuildProcessingImages( InputImagePointer image )
{
  Superclass::BuildProcessingImages( image );

  InputImageListType  lowerImages = this->GetLowerOutlierImages();
  InputImageListType  upperImages = this->GetUpperOutlierImages();

  for( unsigned int i = 0; i < this->GetNumberOfOutlierImagesToRemove();
    i++ )
    {
    InputImagePointer lowerImage = InputImageType::New();
    lowerImage->SetRegions( image->GetLargestPossibleRegion() );
    lowerImage->SetSpacing( image->GetSpacing() );
    lowerImage->SetOrigin( image->GetOrigin() );
    lowerImage->Allocate();
    lowerImage->FillBuffer( NumericTraits<InputPixelType>::max() );

    InputImagePointer upperImage = InputImageType::New();
    upperImage->SetRegions( image->GetLargestPossibleRegion() );
    upperImage->SetSpacing( image->GetSpacing() );
    upperImage->SetOrigin( image->GetOrigin() );
    upperImage->Allocate();
    upperImage->FillBuffer(
      NumericTraits<InputPixelType>::NonpositiveMin() );

    lowerImages.push_back( lowerImage );
    upperImages.push_back( upperImage );
    ::tube::FmtInfoMessage( "Building outlier image %d!", i );
    }

  // Add additional images to the lower outlier list if the median is used
  if( UseMedian() )
    {
    // Number of remaining images to add to lower group to insure that
    // all possible median values are covered
    unsigned int nTotal = this->GetTotalNumberOfImages();
    unsigned int nOutliers = this->GetNumberOfOutlierImagesToRemove();
    unsigned int remaining = ( nTotal / 2 ) + 1 - nOutliers;
    ::tube::FmtInfoMessage( "Images remaining to build: %d!", remaining );

    for( unsigned int i = 0; i < remaining; i++ )
      {
      InputImagePointer lowerImage = InputImageType::New();
      lowerImage->SetRegions( image->GetLargestPossibleRegion() );
      lowerImage->SetSpacing( image->GetSpacing() );
      lowerImage->SetOrigin( image->GetOrigin() );
      lowerImage->Allocate();
      lowerImage->FillBuffer( NumericTraits<InputPixelType>::max() );

      lowerImages.push_back( lowerImage );
      }
    }
  ::tube::FmtInfoMessage( "This is the lower list size: %d",
    lowerImages.size() );

  this->SetLowerOutlierImages( lowerImages );
  this->SetUpperOutlierImages( upperImages );
}

template< class TInputImageType, class TOutputMeanImageType,
  class TOutputSigmaImageType >
void
RobustMeanAndSigmaImageBuilder< TInputImageType, TOutputMeanImageType,
  TOutputSigmaImageType >
::AddImage( InputImagePointer i )
{
  // Add image to running total for the variance calculation
  Superclass::AddImage( i );

  InputImagePointer medianCopy  = GetImageCopy( i );
  InputImagePointer upperCopy   = GetImageCopy( i );

  // Update the median image list
  if( this->UseMedian() )
    {
    // ** Update the lower half of the images
    SetMedianImages( i );
    }
  else
    {
    // Update the lower outlying images
    SetLowerImages( i );
    }

  // Update the upper outlying images
  SetUpperImages( i );
}

template< class TInputImageType, class TOutputMeanImageType,
  class TOutputSigmaImageType >
void
RobustMeanAndSigmaImageBuilder< TInputImageType, TOutputMeanImageType,
  TOutputSigmaImageType >
::SetMedianImages( InputImagePointer i )
{
  // Get image copy, such that input pointer is not changed by fun
  InputImagePointer image = GetImageCopy( i );

  bool isAscending = true;
  AddToUpdateImageList( image, this->GetLowerOutlierImages(), isAscending );

}

template< class TInputImageType, class TOutputMeanImageType,
  class TOutputSigmaImageType >
void
RobustMeanAndSigmaImageBuilder< TInputImageType, TOutputMeanImageType,
  TOutputSigmaImageType >
::SetLowerImages( InputImagePointer i )
{
  bool isAscending = true;
  InputImagePointer image = GetImageCopy( i );
  AddToUpdateImageList( image, this->GetLowerOutlierImages(), isAscending );
}

template< class TInputImageType, class TOutputMeanImageType,
  class TOutputSigmaImageType >
void
RobustMeanAndSigmaImageBuilder< TInputImageType, TOutputMeanImageType,
  TOutputSigmaImageType>
::SetUpperImages( InputImagePointer i )
{
  bool isAscending = false;
  InputImagePointer image = GetImageCopy( i );
  AddToUpdateImageList( image, this->GetUpperOutlierImages(), isAscending );
}

template< class TInputImageType, class TOutputMeanImageType,
  class TOutputSigmaImageType >
typename RobustMeanAndSigmaImageBuilder< TInputImageType,
  TOutputMeanImageType, TOutputSigmaImageType >::InputImageListType&
RobustMeanAndSigmaImageBuilder< TInputImageType, TOutputMeanImageType,
  TOutputSigmaImageType >
::AddToUpdateImageList( InputImagePointer input, InputImageListType& list,
  bool ListIsAscending )
{
  // List is filled from front to back and assumed to be in ascending order
  InputIteratorType  it_input( input, input->GetLargestPossibleRegion() );
  for( int i = 0; i < list.size(); i ++ )
    {
    InputIteratorType  it_list( list[i],
      list[i]->GetLargestPossibleRegion() );

    it_input.GoToBegin();
    it_list.GoToBegin();
    while( !it_list.IsAtEnd() && !it_input.IsAtEnd() )
      {
      // Do not count a pixel if user requests to threshold and the pixel
      // is below the threshold
      if( this->GetThresholdInputImageBelowOn() &&
          it_input.Get() <= this->GetThresholdInputImageBelow() )
        {
        // Do nothing for this pixel ( should not be included in
        // calcualations ( same as itkMeanAndSigmaImageBuilder )
        }
      /* List is ascending & the input value is less than list swap,
       * and vice versa */
      else if( ( ListIsAscending && it_input.Get() < it_list.Get() ) ||
        ( !ListIsAscending && it_input.Get() > it_list.Get() ) )
        {
        // Swap the values and continue
        InputPixelType tmp = it_list.Get();
        it_list.Set( it_input.Get() );
        it_input.Set( tmp );
        }
      ++it_input;
      ++it_list;
      }
    }
  return list;
}

template< class TInputImageType, class TOutputMeanImageType,
  class TOutputSigmaImageType >
typename RobustMeanAndSigmaImageBuilder< TInputImageType,
  TOutputMeanImageType, TOutputSigmaImageType >::InputImagePointer
RobustMeanAndSigmaImageBuilder< TInputImageType, TOutputMeanImageType,
  TOutputSigmaImageType >
::GetImageCopy( InputImagePointer input )
{
  InputImagePointer copy = InputImageType::New();
  copy->SetRegions( input->GetLargestPossibleRegion() );
  copy->SetSpacing( input->GetSpacing() );
  copy->SetOrigin( input->GetOrigin() );
  copy->Allocate();

  InputConstIteratorType it_input( input, input->
    GetLargestPossibleRegion() );
  InputIteratorType it_copy( copy, copy->GetLargestPossibleRegion() );

  while( !it_input.IsAtEnd() )
    {
    it_copy.Set( it_input.Get() );

    ++it_input;
    ++it_copy;
    }
  return copy;
}

template< class TInputImageType, class TOutputMeanImageType,
  class TOutputSigmaImageType >
void
RobustMeanAndSigmaImageBuilder< TInputImageType, TOutputMeanImageType,
  TOutputSigmaImageType >
::FinalizeOutput( void )
{
  if( !( this->GetIsProcessing() ) )
    {
    ::tube::FmtErrorMessage( "Must call Start() before finalizing!" );
    return;
    }

  InputImageListType  lowInputs   = this->GetLowerOutlierImages();
  InputImageListType  highInputs  = this->GetUpperOutlierImages();

  ProcessImagePointer sumImage        = this->GetSumImage();
  ProcessImagePointer sumSquareImage  = this->GetSumSquareImage();
  CountImagePointer   validImages     = this->GetValidCountImage();

  // Output region is defined by the size of the images ( i.e., the largest
  // possible region )
  ProcessIteratorType it_sum( sumImage, sumImage->
    GetLargestPossibleRegion() );
  ProcessIteratorType it_sumSqr( sumSquareImage, sumSquareImage->
    GetLargestPossibleRegion() );
  CountIteratorType   it_valid( validImages, validImages->
    GetLargestPossibleRegion() );

  unsigned int outlierImage = this->GetNumberOfOutlierImagesToRemove();
  for( unsigned int i = 0; i < outlierImage; i++ )
    {
    it_sum.GoToBegin();
    it_sumSqr.GoToBegin();
    it_valid.GoToBegin();

    InputConstIteratorType it_lowInput( lowInputs[i],
      lowInputs[i]->GetLargestPossibleRegion() );
    InputConstIteratorType it_highInput( highInputs[i],
      highInputs[i]->GetLargestPossibleRegion() );

    while( !it_lowInput.IsAtEnd() )
      {
      if( it_valid.Get() > 2*outlierImage )
        {
        ProcessPixelType  sumValue    = it_lowInput.Get() +
          it_highInput.Get();
        ProcessPixelType  sumSqrValue = it_lowInput.Get() *
          it_lowInput.Get() + it_highInput.Get() * it_highInput.Get();
        it_sum.Set( it_sum.Get() - sumValue );
        it_sumSqr.Set( it_sumSqr.Get() - sumSqrValue );
        it_valid.Set( it_valid.Get() - 2 );
        }

      ++it_lowInput;
      ++it_highInput;
      ++it_sum;
      ++it_sumSqr;
      ++it_valid;
      }
    }
  SetSumImage( sumImage );
  SetSumSquareImage( sumSquareImage );
  SetValidCountImage( validImages );

  // Run the finalization using the superclass ( NOTE: Must occur AFTER the
  // subtraction of the outlier images from the summed images )
  Superclass::FinalizeOutput();

  // NOTE: Must occur after the Superclass call to FinalizeOutput() --
  // Otherwise Median will be overwritten by mean
  if( UseMedian() )
    {
    // Build the median image and replace the mean with the median
    SetOutputMeanImage( this->GetMedianImage() );
    }
}

template< class TInputImageType, class TOutputMeanImageType,
  class TOutputSigmaImageType >
typename RobustMeanAndSigmaImageBuilder< TInputImageType,
  TOutputMeanImageType, TOutputSigmaImageType>::OutputMeanImagePointer
RobustMeanAndSigmaImageBuilder< TInputImageType, TOutputMeanImageType,
  TOutputSigmaImageType>
::GetMedianImage( void )
{
  unsigned int totalNumImages = this->GetTotalNumberOfImages();
  InputImageListType lowerImages = this->GetLowerOutlierImages();

  unsigned int numImages = totalNumImages/2;
  typename InputImageType::ConstPointer lastLowerImage =
    lowerImages[numImages];

  // Build output median image
  OutputMeanImagePointer medianImage = OutputMeanImageType::New();
  medianImage->SetRegions( lastLowerImage->GetLargestPossibleRegion() );
  medianImage->SetSpacing( lastLowerImage->GetSpacing() );
  medianImage->SetOrigin( lastLowerImage->GetOrigin() );
  medianImage->Allocate();

  OutputMeanIteratorType it_median( medianImage, medianImage->
    GetLargestPossibleRegion() );
  // Odd number
  if( ( totalNumImages % 2 ) )
    {
    InputConstIteratorType it_middle( lastLowerImage, lastLowerImage->
      GetLargestPossibleRegion() );
    it_middle.GoToBegin();
    it_median.GoToBegin();

    while( !it_middle.IsAtEnd() )
      {
      it_median.Set( it_middle.Get() );
      ++it_middle;
      ++it_median;
      }
    }
  // Even number
  else
    {
    InputConstIteratorType it_middle1( lastLowerImage, lastLowerImage->
      GetLargestPossibleRegion() );
    typename InputImageType::ConstPointer prevLowerImage =
      lowerImages[numImages - 1];
    InputConstIteratorType it_middle2( prevLowerImage,
      prevLowerImage->GetLargestPossibleRegion() );

    it_middle1.GoToBegin();
    it_middle2.GoToBegin();
    it_median.GoToBegin();
    while( !it_middle1.IsAtEnd() && !it_middle2.IsAtEnd() )
      {
      // Average the values
      OutputMeanPixelType median =  ( double( it_middle1.Get() ) +
                                      double( it_middle2.Get() ) ) / 2;
      it_median.Set( median );

      ++it_middle1;
      ++it_middle2;
      ++it_median;
      }
    }
  return medianImage;
}

template< class TInputImageType, class TOutputMeanImageType,
  class TOutputSigmaImageType >
void
RobustMeanAndSigmaImageBuilder< TInputImageType, TOutputMeanImageType,
  TOutputSigmaImageType>
::UpdateOutputImageSize( SizeType inputSize )
{
  if( !( this->GetIsProcessing() ) )
    {
    ::tube::FmtErrorMessage( " Must call AddImage() before updating size" );
    return;
    }

  Superclass::UpdateOutputImageSize( inputSize );

  typedef ResampleImageFilter<InputImageType, InputImageType>
    ResampleInputImageType;

  InputImageListType  lowerList = this->GetLowerOutlierImages();

  typename InputImageListType::const_iterator  it_lower = lowerList.begin();
  while( it_lower != lowerList.end() )
    {
    typename ResampleInputImageType::Pointer filter =
      ResampleInputImageType::New();
    filter->SetInput( *it_lower );
    filter->SetSize( inputSize );
    filter->SetOutputSpacing( ( *it_lower )->GetSpacing() );
    filter->SetOutputOrigin( ( *it_lower )->GetOrigin() );
    filter->Update();

    ++it_lower;
    }

  InputImageListType  upperList = this->GetUpperOutlierImages();

  typename InputImageListType::const_iterator  it_upper = upperList.begin();
  while( it_upper != upperList.end() )
    {
    typename ResampleInputImageType::Pointer filter =
      ResampleInputImageType::New();
    filter->SetInput( *it_upper );
    filter->SetSize( inputSize );
    filter->SetOutputSpacing( ( *it_upper )->GetSpacing() );
    filter->SetOutputOrigin( ( *it_upper )->GetOrigin() );
    filter->Update();

    ++it_upper;
    }
}

#endif // End !defined( __itktubeRobustMeanAndSigmaImageBuilder_hxx )
