
namespace itk
{

namespace tube
{


template<class TInputImage>
MinimizeImageSizeFilter<TInputImage>
::MinimizeImageSizeFilter()
{
  this->SetBufferImage( false );
  m_ThresholdValue = 0;
  this->SetThresholdAbove( false );

  this->SetDefaultPixelValue( 0 );
}



template<class TInputImage>
void
MinimizeImageSizeFilter<TInputImage>
::Update()
{
  // ** Start with complete region ** //
  InputImageConstPointer input = this->GetInput();
  RegionType  region = input->GetLargestPossibleRegion();

  if( TDimension == 3 )
    {
    if( this->GetClipEndIndices() )  // Determine the clipped end of the image
      {
      Get3DCroppedEndRegion( input, region );
      }
    if( this->GetClipStartIndices() )  // Determine the clipped start of the image
      {
      Get3DCroppedStartRegion( input, region );
      }
    }
  else if( TDimension == 2 )
    {
    std::cout << " 2D dimension image cropping is not implemented yet...sorry!" << std::endl;
    // ** (though is easily implemented by copying code for 3D case, removing one imbedded loop and using
    // **    itkLineConstIteratorWithIndex.h instead of itkSlice... ) ** //
    return;
    }
  else
    {
    std::cout << "The dimension size is not compatiable with the filter at this time" << std::endl;
    return;
    }

//std::cout << "New Region: \n" << region << std::endl;

  PointType newOrigin;
  input->TransformIndexToPhysicalPoint( region.GetIndex(), newOrigin );

//std::cout << "This is the final index " << region.GetIndex() << std::endl;
//std::cout << "This is the final origin " << newOrigin << std::endl;

  SizeType size;
  size = region.GetSize();

  // ** Adjust for the buffer (if provided) ** //
  if( this->GetBufferImage() )
    {
    if( this->GetClipEndIndices() )
      {
      for( int i = 0; i < TDimension; i++ )
        {
        size[i] += m_NumberOfBufferPixels[i];
        }
      }
    if( this->GetClipStartIndices() )
      {
      for( int j = 0; j < TDimension; j++ )
        {
        size[j] += m_NumberOfBufferPixels[j];
        // ** Get the physical point for the new origin ** //
        newOrigin[j] -= ( m_NumberOfBufferPixels[j] * input->GetSpacing()[j] );
        }
      }
    }
//std::cout << "New Size: \n" << size << std::endl;

  // ** Resample Image to given region parameters ** //
  typedef ResampleImageFilter<InputImageType,OutputImageType>         ResampleFilterType;
  typename ResampleFilterType::Pointer filter = ResampleFilterType::New();
  filter->SetInput( input );
  filter->SetOutputSpacing( input->GetSpacing() );  // Spacing remains the same

  filter->SetDefaultPixelValue( this->GetDefaultPixelValue() );

  // ** Apply the new size and origin ** //
  filter->SetOutputOrigin( newOrigin );
  filter->SetSize( size );
  filter->Update();

  this->SetOutput( filter->GetOutput() );
}


template<class TInputImage>
void
MinimizeImageSizeFilter<TInputImage>
::Get3DCroppedEndRegion( InputImageConstPointer input, RegionType& region )
{
  InputPixelType  threshold = GetThresholdValue();
  bool            isThresAbove = this->GetThresholdAbove();

  typedef ImageSliceConstIteratorWithIndex< InputImageType>     ImageSliceConstIteratorType;
  SizeType  size = region.GetSize();

//std::cout << "\n\n3DCroppedEndRegion: Size: \n" << size << std::endl;

  // ** Run sweep for back end
  for( unsigned int i = 0; i < TDimension; i++ )
    {
    ImageSliceConstIteratorType  it_input( input, region );
    it_input.SetFirstDirection( i );
    it_input.SetSecondDirection( ( (i+1)%TDimension ) );

    it_input.GoToReverseBegin();

    bool breakPoint = false;
    while( !it_input.IsAtReverseEnd() )
      {
      while( !it_input.IsAtReverseEndOfSlice() )
        {
        while( !it_input.IsAtReverseEndOfLine() )
          {
          if( it_input.Get() > threshold && !isThresAbove ){ breakPoint = true; break; }
          else if( it_input.Get() < threshold && isThresAbove ){ breakPoint = true; break; }

          --it_input;
          }
        if( breakPoint ){ break; }
        it_input.PreviousLine();
        }
      if( breakPoint ){ break; }
      it_input.PreviousSlice();
      }

    // ** This is the index of the first point found in the plane defined by i & i+1 ( starting from end )
    // **   Therefore it is the first point for dimension i+2 ** //
    IndexType index = it_input.GetIndex();

    unsigned int dim = (i+2)%TDimension;  //Dimension that is changed
    if( index[dim] < 0 ) { std::cout << "Index below 0 = " << index[dim] << std::endl; index[dim] = 0; }

    size[dim] = ( index[dim] - region.GetIndex()[dim] ) + 1;
//std::cout << "3DCroppedEndRegion: New index: \n" << index << std::endl;
//std::cout << "3DCroppedEndRegion: New Size: \n" << size << std::endl;

    region.SetSize( size );
    }
}

template<class TInputImage>
void
MinimizeImageSizeFilter<TInputImage>
::Get3DCroppedStartRegion( InputImageConstPointer input, RegionType& region )
{
  InputPixelType  threshold = GetThresholdValue();
  bool            isThresAbove = this->GetThresholdAbove();

  typedef ImageSliceConstIteratorWithIndex< InputImageType>     ImageSliceConstIteratorType;
  IndexType index = region.GetIndex();
  SizeType  size = region.GetSize();

//std::cout << "\n\n3DCroppedStartRegion: Region: \n" << region << std::endl;

  // ** Run sweep for front end of image
  for( unsigned int i = 0; i < TDimension; i++ )
    {
    ImageSliceConstIteratorType  it_input( input, region );
    it_input.SetFirstDirection( i );
    it_input.SetSecondDirection( ( (i+1)%TDimension ) );

    it_input.GoToBegin();

    bool breakPoint = false;
    while( !it_input.IsAtEnd() )
      {
      while( !it_input.IsAtEndOfSlice() )
        {
        while( !it_input.IsAtEndOfLine() )
          {
          if( it_input.Get() > threshold && !isThresAbove ){ breakPoint = true; break; }
          else if( it_input.Get() < threshold && isThresAbove ){ breakPoint = true; break; }
          ++it_input;
          }
        if( breakPoint ){ break; }
        it_input.NextLine();
        }
      if( breakPoint ){ break; }
      it_input.NextSlice();
      }

    // ** This is the index of the first point found in the plane defined by i & i+1 ( starting from beginning )
    // **   Therefore it is the first point for dimension i+2 ** //
    unsigned int dim = (i+2)%TDimension;  //Dimension that is changed
    index[dim] = it_input.GetIndex()[dim];

    size[dim]  -= ( index[dim] - region.GetIndex()[dim] );
    region.SetSize( size );
    region.SetIndex( index );
    }
}

} // End of namespace tube

} // End of itk namespace
