/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeMergeAdjacentImagesFilter_h
#define __itktubeMergeAdjacentImagesFilter_h

// standard includes
#include <vector>

// ITK includes
#include <itkMacro.h>
#include <itkImageSource.h>

namespace itk
{

namespace tube
{

template< typename TImage >
class MergeAdjacentImagesFilter :
  public ImageSource< TImage >
{
public:

  typedef MergeAdjacentImagesFilter                          Self;
  typedef ImageSource< TImage >                              Superclass;
  typedef SmartPointer< Self >                               Pointer;
  typedef SmartPointer< const Self >                         ConstPointer;

  typedef TImage                                             ImageType;
  typedef typename TImage::PixelType                         PixelType;
  typedef std::vector< int >                                 PaddingType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( MergeAdjacentImagesFilter, ImageSource );

  itkStaticConstMacro( ImageDimension, unsigned int, TImage::ImageDimension );

  /** Set input image 1 */
  virtual void SetInput1( const ImageType * image );

  /** Get input image 1 */
  itkGetConstObjectMacro( Input1, ImageType );

  /** Set input image 2 */
  virtual void SetInput2( const ImageType * image );

  /** Get input image 2 */
  itkGetConstObjectMacro( Input2, ImageType );

  /** Set value used for output pixels that dont intersect with input image */
  itkSetMacro( Background, PixelType );

  /** Get value used for output pixels that dont intersect with input image */
  itkGetMacro( Background, PixelType );

  /** Set if zero-valued input pixels should be ignored */
  itkSetMacro( MaskZero, bool );

  /** Get if zero-valued input pixels should be ignored */
  itkGetMacro( MaskZero, bool );

  /** Set number of registration iterations */
  itkSetMacro( MaxIterations, unsigned int );

  /** Get number of registration iterations */
  itkGetMacro( MaxIterations, unsigned int );

  /** Set padding for second image */
  void SetPadding( const PaddingType & padding );

  /** Get padding for second image */
  itkGetConstReferenceMacro( Padding, PaddingType );

  /** Set expected initial misalignment offset */
  itkSetMacro( ExpectedOffset, double );

  /** Get expected initial misalignment offset */
  itkGetMacro( ExpectedOffset, double );

  /** Set expected initial misalignment rotation */
  itkSetMacro( ExpectedRotation, double );

  /** Get expected initial misalignment rotation */
  itkGetMacro( ExpectedRotation, double );

  /** Set portion of pixels to use to compute similarity in registration */
  itkSetMacro( SamplingRatio, double );

  /** Get portion of pixels to use to compute similarity in registration */
  itkGetMacro( SamplingRatio, double );

  /** Set if overlapping pixels should be averaged insted of blending */
  itkSetMacro( BlendUsingAverage, bool );

  /** Get if overlapping pixels should be averaged instead of blending */
  itkGetMacro( BlendUsingAverage, bool );

  /** Set use of experimental method for fast blending */
  itkSetMacro( UseFastBlending, bool );

  /** Get use of experimental method for fast blending */
  itkGetMacro( UseFastBlending, bool );

  /** Set filename to load the transform from */
  void LoadTransform( const std::string & filename );

  /** Set filename to save the transform to */
  void SaveTransform( const std::string & filename );

protected:
  MergeAdjacentImagesFilter( void );
  virtual ~MergeAdjacentImagesFilter( void ) {}

  virtual void GenerateData() override;

  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:
  // Purposely not implemented
  MergeAdjacentImagesFilter( const Self & );
  void operator = ( const Self & );

  typename TImage::PixelType             m_Background;
  bool                                   m_MaskZero;
  unsigned int                           m_MaxIterations;
  PaddingType                            m_Padding;
  double                                 m_ExpectedOffset;
  double                                 m_ExpectedRotation;
  double                                 m_SamplingRatio;
  bool                                   m_BlendUsingAverage;
  bool                                   m_UseFastBlending;
  std::string                            m_InitialTransformFile;
  std::string                            m_OutputTransformFile;

  typename ImageType::ConstPointer       m_Input1;
  typename ImageType::ConstPointer       m_Input2;

};  // End class MergeAdjacentImagesFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeMergeAdjacentImagesFilter.hxx"
#endif

#endif
