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

#ifndef __itktubeSpatialObjectToImageFilter_h
#define __itktubeSpatialObjectToImageFilter_h

#include <itkSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkTubeSpatialObject.h>
#include <itkTubeSpatialObjectPoint.h>

namespace itk
{

namespace tube
{

/** \class TubeSpatialObjectToImageFilter
 * \brief This filter creates a binary image with 1 representing the
 * vessel existence in that voxels and 0 not.
 * Also, forms the same image, but with the radius value in place of the 1.
 */

template< unsigned int ObjectDimension, class TOutputImage,
  class TRadiusImage = Image< float, TOutputImage::ImageDimension >,
  class TTangentImage = Image< Vector< float, TOutputImage::ImageDimension >,
                                       TOutputImage::ImageDimension > >
class TubeSpatialObjectToImageFilter
  : public SpatialObjectToImageFilter< SpatialObject< ObjectDimension >,
                                       TOutputImage >
{
public:

  /** Standard class typedefs. */
  typedef TubeSpatialObjectToImageFilter                 Self;
  typedef SpatialObjectToImageFilter< SpatialObject<ObjectDimension>,
                                      TOutputImage>      SuperClass;
  typedef SmartPointer< Self >                           Pointer;
  typedef SmartPointer< const Self >                     ConstPointer;

  /** Tube class typedef */
  typedef TOutputImage                                   OutputImageType;
  typedef SpatialObject<ObjectDimension>                 SpatialObjectType;
  typedef typename SpatialObjectType::ChildrenListType   ChildrenListType;
  typedef TubeSpatialObject<ObjectDimension>             TubeType;
  typedef typename TOutputImage::SizeType                SizeType;

  typedef TRadiusImage                                   RadiusImage;
  typedef typename TRadiusImage::Pointer                 RadiusImagePointer;
  typedef typename TRadiusImage::PixelType               RadiusPixelType;

  typedef TTangentImage                                  TangentImage;
  typedef typename TTangentImage::Pointer                TangentImagePointer;
  typedef typename TTangentImage::PixelType              TangentPixelType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( TubeSpatialObjectToImageFilter,
                SpatialObjectToImageFilter );

  /** Set if the tube should be full inside */
  itkSetMacro( UseRadius, bool );
  itkGetMacro( UseRadius, bool );

  /** Set if the filter should build a radius image in parallel */
  itkSetMacro( BuildRadiusImage, bool );
  itkGetMacro( BuildRadiusImage, bool );

  /** Set if the filter should build a tangent image in parallel */
  itkSetMacro( BuildTangentImage, bool );
  itkGetMacro( BuildTangentImage, bool );

  /**Image Pointer Definition and GetRadius return method */
  RadiusImagePointer GetRadiusImage( void );

  /**Image Pointer Definition and GetRadius return method */
  TangentImagePointer GetTangentImage( void );

  /** Set the FallOff value */
  itkSetMacro( FallOff, double );
  itkGetMacro( FallOff, double );

  /** Set if the value of tubes that are crossing should accumulate
   *  their values to produce the image */
  itkSetMacro( Cumulative, bool );
  itkGetMacro( Cumulative, bool );

protected:

  TubeSpatialObjectToImageFilter( void );
  ~TubeSpatialObjectToImageFilter( void );

  /** Create the output images and fill it */
  void GenerateData( void );

  void PrintSelf(std::ostream& os, Indent indent) const
    {
    SuperClass::PrintSelf(os,indent);
    os << indent << "m_UseRadius: " << m_UseRadius << std::endl;
    os << indent << "m_FallOff: " << m_FallOff << std::endl;
    os << indent << "m_Cumulative: " << m_Cumulative << std::endl;
    }

private:

  bool        m_BuildRadiusImage;
  bool        m_BuildTangentImage;
  bool        m_UseRadius;
  double      m_FallOff;
  bool        m_Cumulative;

  typename RadiusImage::Pointer     m_RadiusImage;
  typename TangentImage::Pointer    m_TangentImage;

}; // End class TubeSpatialObjectToImageFilter

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSpatialObjectToImageFilter.hxx"
#endif

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeSpatialObjectToImageFilter_h)
