/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __itktubeConvertShrunkenSeedImageToListFilter_h
#define __itktubeConvertShrunkenSeedImageToListFilter_h

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
  /** Standard class typedefs. */
  typedef ConvertShrunkenSeedImageToListFilter             Self;
  typedef ProcessObject                                    Superclass;
  typedef SmartPointer< Self >                             Pointer;
  typedef SmartPointer< const Self >                       ConstPointer;

  typedef TImage                                           ImageType;
  typedef typename ImageType::Pointer                      ImagePointer;
  typedef typename ImageType::PixelType                    PixelType;
  typedef typename ImageType::IndexType                    IndexType;
  typedef ImageFileReader< ImageType >                     ReaderType;
  typedef ImageRegionIterator< ImageType >                 IteratorType;

  typedef TPointsImage                                     PointsImageType;
  typedef typename PointsImageType::Pointer                PointsImagePointer;
  typedef typename PointsImageType::PixelType              PointsPixelType;
  typedef typename PointsImageType::IndexType              PointsIndexType;
  typedef ImageFileReader< PointsImageType >               PointsReaderType;
  typedef ImageRegionIterator< PointsImageType >           PointsIteratorType;

  typedef vnl_matrix< PixelType >                          VnlMatrixType;
  typedef SimpleDataObjectDecorator< VnlMatrixType >       OutputType;

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
  virtual void GenerateData() ITK_OVERRIDE;
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;
  void VerifyPreconditions();

private:
  /** itkConvertShrunkenSeedImageToListFilter parameters **/
  ConvertShrunkenSeedImageToListFilter( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const Superclass::DataObjectIdentifierType &,
    itk::DataObject * ) {};

  VnlMatrixType                        m_VnlOutput;
  double                               m_Threshold;

}; // End class ConvertShrunkenSeedImageToListFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeConvertShrunkenSeedImageToListFilter.hxx"
#endif

#endif // End !defined(_itktubeConvertShrunkenSeedImageToListFilter_h)
