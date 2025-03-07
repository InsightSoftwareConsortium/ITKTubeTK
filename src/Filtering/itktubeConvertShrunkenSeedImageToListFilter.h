/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef itktubeConvertShrunkenSeedImageToListFilter_h
#define itktubeConvertShrunkenSeedImageToListFilter_h

#include <itkProcessObject.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>
#include <itkSimpleDataObjectDecorator.h>

namespace itk
{
namespace tube
{
/** \class ConvertShrunkenSeedImageToList
 */

template< class TImage, class TPointsImage >
class ConvertShrunkenSeedImageToListFilter
  : public ProcessObject
{
public:
  /** Standard class type alias. */
  using Self = ConvertShrunkenSeedImageToListFilter;
  using Superclass = ProcessObject;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  using ImageType = TImage;
  using ImagePointer = typename ImageType::Pointer;
  using PixelType = typename ImageType::PixelType;
  using IndexType = typename ImageType::IndexType;
  using ReaderType = ImageFileReader< ImageType >;
  using IteratorType = ImageRegionIterator< ImageType >;

  using PointsImageType = TPointsImage;
  using PointsImagePointer = typename PointsImageType::Pointer;
  using PointsPixelType = typename PointsImageType::PixelType;
  using PointsIndexType = typename PointsImageType::IndexType;
  using PointsReaderType = ImageFileReader< PointsImageType >;
  using PointsIteratorType = ImageRegionIterator< PointsImageType >;

  using VnlMatrixType = vnl_matrix< PixelType >;
  using OutputType = SimpleDataObjectDecorator< VnlMatrixType >;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ConvertShrunkenSeedImageToListFilter, ProcessObject );

  /** ImageDimension constants */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TImage::ImageDimension );

  /** Method to set/get the image */
  void SetInput( const ImageType* image );
  const ImageType* GetInput( void ) const;

  /** Method to set/get the scale image */
  void SetScaleImage( const ImageType* image );
  const ImageType* GetScaleImage( void ) const;

  /** Method to set/get the points image */
  void SetPointsImage( const PointsImageType* image );
  const PointsImageType* GetPointsImage( void ) const;

  SimpleDataObjectDecorator< vnl_matrix <typename TImage::PixelType> >*
    GetOutput();

  itkGetMacro( Threshold, double );
  itkSetMacro( Threshold, double );

protected:
  ConvertShrunkenSeedImageToListFilter( void );
  ~ConvertShrunkenSeedImageToListFilter( void ) {};

  virtual void GenerateData() override;

  void PrintSelf( std::ostream & os, itk::Indent indent ) const override;

  void VerifyPreconditions();

private:
  /** itkConvertShrunkenSeedImageToListFilter parameters **/
  ConvertShrunkenSeedImageToListFilter( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const typename Superclass::DataObjectIdentifierType &,
    itk::DataObject * ) override {};

  VnlMatrixType                        m_VnlOutput;
  double                               m_Threshold;

}; // End class ConvertShrunkenSeedImageToListFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeConvertShrunkenSeedImageToListFilter.hxx"
#endif

#endif // End !defined(_itktubeConvertShrunkenSeedImageToListFilter_h)
