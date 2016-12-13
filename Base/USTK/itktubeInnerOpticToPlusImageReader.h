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

#ifndef __itktubeInnerOpticToPlusImageReader_h
#define __itktubeInnerOpticToPlusImageReader_h

#include <itkImageSource.h>
#include <itkRGBPixel.h>

#include "SyncRecordManager.h"

namespace itk
{

namespace tube
{

/** \class InnerOpticToPlusImageReader
 *
 * \brief Convert screen captured, tracked InnerOptic frames to 3D Plus Image.
 *
 * Convert the tracking metadata and recorded frames in the InnerOptic
 * extended .ppm files into an itk::Image in the Plus ultrasound library
 * format.
 *
 * To extract only a subset of the images referenced in the InnerOptic
 * metadata file, use SetStartIndex, SetEndIndex, and Set IncrementIndex.
 *
 */
class InnerOpticToPlusImageReader
  : public ImageSource< Image< RGBPixel< unsigned char >, 3 > >
{
public:
  static const unsigned int ImageDimension = 3;

  typedef unsigned char                      PixelComponentType;
  typedef RGBPixel< PixelComponentType >     PixelType;
  typedef Image< PixelType, ImageDimension > OutputImageType;

  /** Standard class typedefs. */
  typedef InnerOpticToPlusImageReader    Self;
  typedef ImageSource< OutputImageType > Superclass;
  typedef SmartPointer< Self >           Pointer;
  typedef SmartPointer< const Self >     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( InnerOpticToPlusImageReader, ImageSource );

  /** Set the name of the InnerOptic metadata file for this series. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );

  /** Set the start index of the frames to extract. */
  itkSetMacro( StartIndex, SizeValueType );
  itkGetConstMacro( StartIndex, SizeValueType );

  /** Set the end index of the frames to extract. */
  itkSetMacro( EndIndex, SizeValueType );
  itkGetConstMacro( EndIndex, SizeValueType );

  /** Set in the increment of the frames to extract.  The default value is 1. */
  itkSetMacro( IncrementIndex, SizeValueType );
  itkGetConstMacro( IncrementIndex, SizeValueType );


protected:
  InnerOpticToPlusImageReader( void );
  virtual ~InnerOpticToPlusImageReader( void );

  typedef Superclass::OutputImageRegionType OutputImageRegionType;

  virtual void GenerateOutputInformation( void );

  virtual void GenerateData( void );

  /** Expand to the largest value X and Y -- leave Z unchanged. */
  virtual void EnlargeOutputRequestedRegion( DataObject * output );

private:
  InnerOpticToPlusImageReader( const Self & ); // purposely not implemented
  void operator=( const Self & ); // purposely not implemented

  std::string m_FileName;

  SyncRecordManager * m_SyncRecordManager;

  SizeValueType m_StartIndex;
  SizeValueType m_EndIndex;
  SizeValueType m_IncrementIndex;

}; // End class InnerOpticToPlusImageReader

} // End namespace tubetk

} // End namespace itk

#endif // End !defined( __itktubeInnerOpticToPlusImageReader_h )
