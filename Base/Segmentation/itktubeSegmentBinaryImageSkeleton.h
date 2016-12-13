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

#ifndef __itktubeSegmentBinaryImageSkeleton_h
#define __itktubeSegmentBinaryImageSkeleton_h

// ITK includes
#include <itkBinaryThinningImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>

namespace itk
{

namespace tube
{

/** \class SegmentBinaryImageSkeleton
 * \brief Computes skeleton of a binary image.
 * The output skeleton can be dilated if a radius greater than zero is
 * provided
 */

template< class TPixel, unsigned int VDimension >
class SegmentBinaryImageSkeleton
  : public itk::BinaryThinningImageFilter
      <itk::Image< TPixel, VDimension >, itk::Image< TPixel, VDimension > >
{
public:
  /** custom typedefs */
  typedef itk::Image< TPixel, VDimension > ImageType;
  typedef itk::BinaryBallStructuringElement< TPixel, VDimension>
                                           SEType;
  typedef itk::BinaryDilateImageFilter< ImageType, ImageType, SEType >
                                           DilateType;

  /** Standard class typedefs. */
  typedef SegmentBinaryImageSkeleton                             Self;
  typedef itk::BinaryThinningImageFilter< ImageType, ImageType > Superclass;
  typedef SmartPointer< Self >                                   Pointer;
  typedef SmartPointer< const Self >                             ConstPointer;


  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( SegmentBinaryImageSkeleton, BinaryThinningImageFilter );

  /** Set/Get radius for post-dilatation */
  itkSetMacro( Radius, unsigned int );
  itkGetMacro( Radius, unsigned int );

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro( InputConvertibleToOutputCheck,
    ( Concept::Convertible< TPixel, unsigned char > ) );
  // End concept checking
#endif

protected:

  SegmentBinaryImageSkeleton( void );
  ~SegmentBinaryImageSkeleton( void ) {};

  void PrintSelf( std::ostream& os, Indent indent ) const;

  void GenerateData( void );

private:

  unsigned int m_Radius;

}; // End class SegmentBinaryImageSkeleton

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSegmentBinaryImageSkeleton.hxx"
#endif

#endif // End !defined( __itktubeSegmentBinaryImageSkeleton_h )
