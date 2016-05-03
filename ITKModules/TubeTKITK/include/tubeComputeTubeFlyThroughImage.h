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
#ifndef __tubeComputeTubeFlyThroughImage_h
#define __tubeComputeTubeFlyThroughImage_h

// ITK includes
#include <itkProcessObject.h>
#include <itkGroupSpatialObject.h>

// TubeTK includes
#include "itktubeComputeTubeFlyThroughImageFilter.h"

namespace tube
{
/** \class ComputeTubeFlyThroughImage
 *
 *  \ingroup TubeTKITK
 */

template< class TPixel, unsigned int Dimension >
class ComputeTubeFlyThroughImage:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ComputeTubeFlyThroughImage                 Self;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  typedef TubeSpatialObject< Dimension >             TubesType;
  typedef Image< TPixel, Dimension >                 InputImageType;
  typedef InputImageType                             OutputImageType;
  typedef Image< unsigned char, Dimension >          OutputMaskType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ComputeTubeFlyThroughImage, Object);

  /** Set tube id for which the fly through image is to be generated */
  void SetTubeId( unsigned long tubeId );
  unsigned long GetTubeId();

  /* Set input image from which the tubes were extracted/segmented */
  void SetInputImage( typename InputImageType::Pointer inputImage );

  /* Set tubes */
  void SetInput( typename TubesType::Pointer pTubes );

  /* Generates tube fly through image and mask */
  void Update();

  /* Get the generated tube fly through image */
  typename OutputImageType::Pointer GetOutput();

  /* Get the generated tube fly through image */
  typename OutputMaskType::Pointer GetOutputMask();

protected:
  ComputeTubeFlyThroughImage( void );
  ~ComputeTubeFlyThroughImage() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itkComputeTubeFlyThroughImageFilter parameters **/
  ComputeTubeFlyThroughImage(const Self &);
  void operator=(const Self &);

  typedef itk::tube::TubeSpatialObjectToImageFilter< Dimension,
    OutputImageType > TubesToImageFilterType;

  typename TubesToImageFilterType::Pointer m_TubesToImageFilter;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeComputeTubeFlyThroughImage.hxx"
#endif

#endif // End !defined( __tubeComputeTubeFlyThroughImage_h )
