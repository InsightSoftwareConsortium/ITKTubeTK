/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/
#ifndef __tubeCropTubes_h
#define __tubeCropTubes_h

#include "itktubeCropTubesFilter.h"
#include "itkObject.h"
#include "tubeWrappingMacros.h"

namespace tube
{
/** \class ComputeBinaryImageSimilarityMetrics
 *
 *  \ingroup TubeTKITK
 */

template< unsigned int VDimension >
class CropTubes:
  public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef CropTubes                                Self;
  typedef itk::SmartPointer< Self >                Pointer;
  typedef itk::SmartPointer< const Self >          ConstPointer;
  typedef itk::tube::CropTubesFilter< VDimension > FilterType;

  typedef typename FilterType::TubeGroupType       TubeGroupType;
  typedef typename TubeGroupType::Pointer          TubeGroupPointer;
  typedef typename FilterType::ImageType           ImageType;
  typedef typename FilterType::PointType           PointType;
  typedef typename FilterType::VectorType          VectorType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( CropTubes, Object );

 /* Set input tubes */
  tubeWrapSetMacro( Input, TubeGroupPointer, Filter );
  tubeWrapGetConstObjectMacro( Input, TubeGroupType, Filter );

  /** Set/Get bounding box corner */
  tubeWrapSetMacro( BoxPosition, PointType, Filter );
  tubeWrapGetMacro( BoxPosition, PointType, Filter );

  /** Set/Get bounding box size */
  tubeWrapSetMacro( BoxSize, VectorType, Filter );
  tubeWrapGetMacro( BoxSize, VectorType, Filter );

  /** Set/Get mask image */
  tubeWrapSetObjectMacro( MaskImage, ImageType, Filter );
  tubeWrapGetObjectMacro( MaskImage, ImageType, Filter );

  /** Set/Get Use mask image to crop tubes */
  tubeWrapSetMacro( UseMaskImage, bool, Filter );
  tubeWrapGetMacro( UseMaskImage, bool, Filter );

  /** Set/Get bool to crop tubes or not */
  tubeWrapSetMacro( CropTubes, bool, Filter );
  tubeWrapGetMacro( CropTubes, bool, Filter );

  /** Run Crop Tubes application */
  tubeWrapUpdateMacro( Filter );

  /* Get the crop tubes output */
  tubeWrapGetMacro( Output, TubeGroupPointer, Filter );

protected:
  CropTubes( void );
  ~CropTubes() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itktubeCropTubesFilter parameters **/
  CropTubes( const Self & );
  void operator=( const Self & );

  typename FilterType::Pointer m_Filter;

};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeCropTubes.hxx"
#endif

#endif // End !defined( __tubeCropTubes_h )
