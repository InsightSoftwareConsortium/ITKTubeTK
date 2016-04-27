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
#ifndef __tubeConvertTubesToImage_h
#define __tubeConvertTubesToImage_h

#include <itkObject.h>
#include <itkSpatialObjectReader.h>

#include "itktubeSpatialObjectToImageFilter.h"

namespace tube
{
/** \class ConvertTubesToImage
 *
 *  \ingroup TubeTKITK
 */

template< class TPixel, unsigned int Dimension >
class ConvertTubesToImage:
  public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef ConvertTubesToImage                        Self;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  typedef TPixel                                     TemplatePixelType;
  typedef unsigned char                              OutputPixelType;

  typedef itk::GroupSpatialObject< Dimension >       TubesType;
  typedef itk::Image< TemplatePixelType >            TemplateImageType;
  typedef itk::Image< OutputPixelType >              OutputImageType;


  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ConvertTubesToImage, Object);

  /** Set if the tube should be full inside */
  itkSetMacro( UseRadius, bool );
  itkGetMacro( UseRadius, bool );
  itkBooleanMacro( UseRadius );

  /* Set template image */
  void SetTemplateImage( TemplateImageType::Pointer pTemplateImage );

  /* Set input tubes */
  void SetInput( TubesType::Pointer pTubes );

  void Update();

  typename OutputImageType::Pointer GetOutput();

protected:
  ConvertTubesToImage( void );
  ~ConvertTubesToImage() {}
  void PrintSelf(std::ostream & os, itk::Indent indent) const;

private:
  /** itkConvertTubesToImageFilter parameters **/
  ConvertTubesToImage(const Self &);
  void operator=(const Self &);

  typedef itk::tube::SpatialObjectToImageFilter< Dimension, OutputImageType >
    TubesToImageFilterType;

  typename TubesToImageFilterType::Pointer m_TubesToImageFilter;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeConvertTubesToImage.hxx"
#endif

#endif // End !defined( __tubeConvertTubesToImage_h )
