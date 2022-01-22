/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#ifndef __itktubeTubeSpatialObjectToImageFilter_h
#define __itktubeTubeSpatialObjectToImageFilter_h

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
  typedef SpatialObjectToImageFilter< SpatialObject< ObjectDimension >,
                                      TOutputImage>      SuperClass;
  typedef SmartPointer< Self >                           Pointer;
  typedef SmartPointer< const Self >                     ConstPointer;

  /** Tube class typedef */
  typedef TOutputImage                                   OutputImageType;

  typedef SpatialObject<ObjectDimension>                 SpatialObjectType;
  typedef typename SpatialObjectType::ChildrenListType   ChildrenListType;
  typedef TubeSpatialObject<ObjectDimension>             TubeType;

  typedef TRadiusImage                                   RadiusImage;
  typedef typename TRadiusImage::Pointer                 RadiusImagePointer;
  typedef typename TRadiusImage::PixelType               RadiusPixelType;

  typedef TTangentImage                                  TangentImage;
  typedef typename TTangentImage::Pointer                TangentImagePointer;
  typedef typename TTangentImage::PixelType              TangentPixelType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
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

  /** Set if the value of tubes that are crossing should accumulate
   *  their values to produce the image */
  itkSetMacro( Cumulative, bool );
  itkGetMacro( Cumulative, bool );

  itkSetMacro( ColorByTubeId, bool );
  itkSetMacro( ColorByPointId, bool );
  itkSetMacro( ColorByRadius, bool );
  itkSetMacro( ColorByRidgeness, bool );
  itkSetMacro( ColorByMedialness, bool );
  itkSetMacro( ColorByBranchness, bool );
  itkSetMacro( ColorByCurvature, bool );
  itkSetMacro( ColorByLevelness, bool );
  itkSetMacro( ColorByRoundness, bool );
  itkSetMacro( ColorByIntensity, bool );

  itkGetMacro( ColorByTubeId, bool );
  itkGetMacro( ColorByPointId, bool );
  itkGetMacro( ColorByRadius, bool );
  itkGetMacro( ColorByRidgeness, bool );
  itkGetMacro( ColorByMedialness, bool );
  itkGetMacro( ColorByBranchness, bool );
  itkGetMacro( ColorByCurvature, bool );
  itkGetMacro( ColorByLevelness, bool );
  itkGetMacro( ColorByRoundness, bool );
  itkGetMacro( ColorByIntensity, bool );

protected:

  TubeSpatialObjectToImageFilter( void );
  ~TubeSpatialObjectToImageFilter( void );

  /** Create the output images and fill it */
  void GenerateData( void ) override;

  void PrintSelf( std::ostream& os, Indent indent ) const override
    {
    SuperClass::PrintSelf( os, indent );
    os << indent << "m_BuildRadiusImage: " << m_UseRadius << std::endl;
    os << indent << "m_BuildTangentImage: " << m_UseRadius << std::endl;
    os << indent << "m_UseRadius: " << m_UseRadius << std::endl;
    os << indent << "m_Cumulative: " << m_Cumulative << std::endl;
    os << indent << "m_ColorByTubeId: " << m_ColorByTubeId << std::endl;
    os << indent << "m_ColorByPointId: " << m_ColorByPointId << std::endl;
    os << indent << "m_ColorByRadius: " << m_ColorByRadius << std::endl;
    os << indent << "m_ColorByRidgeness: " << m_ColorByRidgeness << std::endl;
    os << indent << "m_ColorByMedialness: " << m_ColorByMedialness << std::endl;
    os << indent << "m_ColorByBranchness: " << m_ColorByBranchness << std::endl;
    os << indent << "m_ColorByCurvature: " << m_ColorByCurvature << std::endl;
    os << indent << "m_ColorByLevelness: " << m_ColorByLevelness << std::endl;
    os << indent << "m_ColorByRoundness: " << m_ColorByRoundness << std::endl;
    os << indent << "m_ColorByIntensity: " << m_ColorByIntensity << std::endl;
    }

private:

  bool        m_BuildRadiusImage;
  bool        m_BuildTangentImage;
  bool        m_UseRadius;
  bool        m_Cumulative;

  bool        m_ColorByTubeId;
  bool        m_ColorByPointId;
  bool        m_ColorByRadius;
  bool        m_ColorByRidgeness;
  bool        m_ColorByMedialness;
  bool        m_ColorByBranchness;
  bool        m_ColorByCurvature;
  bool        m_ColorByLevelness;
  bool        m_ColorByRoundness;
  bool        m_ColorByIntensity;

  typename RadiusImage::Pointer     m_RadiusImage;
  typename TangentImage::Pointer    m_TangentImage;

}; // End class TubeSpatialObjectToImageFilter

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTubeSpatialObjectToImageFilter.hxx"
#endif

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeTubeSpatialObjectToImageFilter_h )
