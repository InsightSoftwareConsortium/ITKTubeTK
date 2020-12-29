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

#ifndef __itktubeCropTubesFilter_h
#define __itktubeCropTubesFilter_h

#include "itkGroupSpatialObject.h"
#include "itkImage.h"
#include "itktubeSpatialObjectToSpatialObjectFilter.h"
#include "itkTubeSpatialObject.h"
#include "itkTubeSpatialObjectPoint.h"
namespace itk
{
namespace tube
{

/** \class CropTubesFilter
 * \brief Crop tubes spatial object based on a volume map or a bounding box.
 *
 */
template< unsigned int VDimension >
class CropTubesFilter
  : public SpatialObjectToSpatialObjectFilter<
    GroupSpatialObject< VDimension >, GroupSpatialObject< VDimension > >
{
public:
  /** Standard class typedefs. */
  typedef GroupSpatialObject< VDimension >      TubeGroupType;

  typedef CropTubesFilter                            Self;
  typedef SpatialObjectToSpatialObjectFilter< TubeGroupType, TubeGroupType >
                                                     Superclass;
  typedef SmartPointer< Self >                       Pointer;
  typedef SmartPointer< const Self >                 ConstPointer;
  typedef TubeSpatialObject< VDimension >      TubeType;
  typedef TubeSpatialObjectPoint< VDimension > TubePointType;

  typedef double                                     PixelType;
  typedef itk::Image< PixelType, VDimension >        ImageType;
  typedef itk::Vector< PixelType, VDimension >       VectorType;
  typedef itk::Point< double, VDimension >           PointType;

  /** Run-time type information ( and related methods ).   */
  itkTypeMacro( CropTubesFilter,
                SpatialObjectToSpatialObjectFilter );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Set/Get bounding box corner */
  itkSetMacro( BoxPositionInWorldSpace, PointType );
  itkGetMacro( BoxPositionInWorldSpace, PointType );

  /** Set/Get bounding box size */
  itkSetMacro( BoxSizeInWorldSpace, VectorType );
  itkGetMacro( BoxSizeInWorldSpace, VectorType );

  /** Set/Get mask image */
  itkSetObjectMacro( MaskImage, ImageType );
  itkGetModifiableObjectMacro( MaskImage, ImageType );

  /** Set/Get Use mask image to crop tubes */
  itkSetMacro( UseMaskImage, bool );
  itkGetMacro( UseMaskImage, bool );

  /** Set/Get bool to crop tubes or not */
  itkSetMacro( CropTubes, bool );
  itkGetMacro( CropTubes, bool );

protected:
  CropTubesFilter( void );
  virtual ~CropTubesFilter( void );

  virtual void GenerateData( void ) override;

  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:
  // purposely not implemented
  CropTubesFilter( const Self & );
  // purposely not implemented
  void operator=( const Self & );

  PointType                   m_BoxPositionInWorldSpace;
  VectorType                  m_BoxSizeInWorldSpace;
  typename ImageType::Pointer m_MaskImage;
  bool                        m_CropTubes;
  bool                        m_UseMaskImage;

  bool IsInsideInWorldSpace( itk::Point< double, VDimension > pointPos,
    double tubeRadius,
    itk::Point< double, VDimension > boxPos,
    itk::Vector< double, VDimension > boxSize,
    std::vector<  typename
    itk::TubeSpatialObjectPoint< VDimension >::
    CovariantVectorType > normalList );

}; // End class CropTubesFilter

} // End namespace tube
} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeCropTubesFilter.hxx"
#endif

#endif // End !defined( __itktubeCropTubesFilter_h )
