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
#ifndef __tubeResampleImage_h
#define __tubeResampleImage_h

#include "itktubeResampleImageFilter.h"
#include "tubeWrappingMacros.h"
#include "itkObject.h"

namespace tube
{
/** \class ResampleImage
 *
 *  \ingroup TubeTKITK
 */

template< class TPixel, unsigned int VDimension >
class ResampleImage:
  public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef ResampleImage                   Self;
  typedef itk::SmartPointer< Self >       Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  typedef itk::tube::ResampleImageFilter
    < TPixel, VDimension >           FilterType;

  typedef typename FilterType::ImageType      ImageType;
  typedef typename ImageType::ConstPointer    ConstImagePointer;
  typedef typename ImageType::Pointer         ImagePointer;
  typedef typename FilterType::TransformType  TransformType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ResampleImage, Object );

  /* Set input image */
  tubeWrapSetConstObjectMacro( Input, ImageType, Filter );

  /** Set/Get input Match Image */
  tubeWrapSetObjectMacro( MatchImage, ImageType, Filter );
  tubeWrapGetObjectMacro( MatchImage, ImageType, Filter );

  /** Set/Get whether Output is isotropic or not */
  tubeWrapSetMacro( MakeIsotropic, bool, Filter );
  tubeWrapGetMacro( MakeIsotropic, bool, Filter );

  /** Set/Get whether Output is high resolution isotropic or not */
  tubeWrapSetMacro( MakeHighResIso, bool, Filter );
  tubeWrapGetMacro( MakeHighResIso, bool, Filter );

  /** Set/Get interpolator */
  tubeWrapSetMacro( Interpolator, std::string, Filter );
  tubeWrapGetMacro( Interpolator, std::string, Filter );

  /** Set/Get whether to load transform or not */
  tubeWrapSetMacro( LoadTransform, bool, Filter );
  tubeWrapGetMacro( LoadTransform, bool, Filter );

  tubeWrapGetObjectMacro( Output, ImageType, Filter );

  /** Set Output Transform */
  void SetTransform( TransformType* t );

  /** Set Output Spacing */
  void SetSpacing( std::vector<double> s );

  /** Set Output Origin */
  void SetOrigin( std::vector<double> o );

  /** Set Output Index */
  void SetIndex( std::vector<int> i );

  /** Set Output Resample Factor */
  void SetResampleFactor( std::vector<double> rf );

  /* Runs tubes to image conversion */
  tubeWrapUpdateMacro( Filter );

protected:
  ResampleImage( void );
  ~ResampleImage() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itktubeResampleImageFilter parameters **/
  ResampleImage( const Self & );
  void operator=( const Self & );

  typename FilterType::Pointer  m_Filter;
};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeResampleImage.hxx"
#endif

#endif // End !defined( __tubeResampleImage_h )
