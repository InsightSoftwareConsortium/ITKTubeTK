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

#ifndef __itktubeComputeTubeFlyThroughImageFilter_h
#define __itktubeComputeTubeFlyThroughImageFilter_h

#include <itkGroupSpatialObject.h>
#include <itkSpatialObjectToImageFilter.h>
#include <itkVesselTubeSpatialObject.h>

namespace itk
{

namespace tube
{

/** \class ComputeTubeFlyThroughImageFilter
 * \brief This filter computes the fly through image and mask
 * for a specified tube
 */

template< class TPixel, unsigned int Dimension >
class ComputeTubeFlyThroughImageFilter
  : public SpatialObjectToImageFilter< GroupSpatialObject< Dimension >,
                                       Image< TPixel, Dimension > >
{
public:

  /** Tube class typedef */
  typedef GroupSpatialObject< Dimension >                TubeGroupType;
  typedef VesselTubeSpatialObject< Dimension >           TubeType;
  typedef Image< TPixel, Dimension >                     InputImageType;
  typedef InputImageType                                 OutputImageType;
  typedef Image< unsigned char, Dimension >              OutputMaskType;

  /** Standard class typedefs. */
  typedef ComputeTubeFlyThroughImageFilter               Self;
  typedef SpatialObjectToImageFilter< TubeGroupType,
    InputImageType >                                     SuperClass;
  typedef SmartPointer< Self >                           Pointer;
  typedef SmartPointer< const Self >                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ComputeTubeFlyThroughImageFilter,
                SpatialObjectToImageFilter );

  /** Set/Get tube id for which the fly through image is to be generated */
  itkSetMacro( TubeId, unsigned long );
  itkGetMacro( TubeId, unsigned long );

  /** Set/Get input image from which the tubes were extracted/segmented */
  itkSetConstObjectMacro( InputImage, InputImageType );
  itkGetConstObjectMacro( InputImage, InputImageType );

  /** Get output tube mask image */
  itkGetObjectMacro( OutputMask, OutputMaskType );

protected:

  ComputeTubeFlyThroughImageFilter( void );
  ~ComputeTubeFlyThroughImageFilter( void ) {};

  /** Creates the tube fly through image */
  void GenerateData( void );

  void PrintSelf( std::ostream& os, Indent indent ) const;

private:

  unsigned long                               m_TubeId;
  typename InputImageType::ConstPointer       m_InputImage;
  typename OutputMaskType::Pointer            m_OutputMask;

}; // End class ComputeTubeFlyThroughImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeComputeTubeFlyThroughImageFilter.hxx"
#endif

#endif // End !defined( __itktubeComputeTubeFlyThroughImageFilter_h )
