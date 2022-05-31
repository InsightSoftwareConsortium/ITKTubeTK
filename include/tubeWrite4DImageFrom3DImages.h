/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __tubeWrite4DImageFrom3DImages_h
#define __tubeWrite4DImageFrom3DImages_h

// ITK Includes
#include "itkProcessObject.h"
#include "itkGroupSpatialObject.h"
#include "itkTubeSpatialObject.h"

// TubeTK Includes
#include "tubeWrappingMacros.h"
#include <itkImage.h>

namespace tube
{
/** \class Write4DImageFrom3DImages
 *
 *  \ingroup TubeTK
 */

template< class InputImageT >
class Write4DImageFrom3DImages:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef Write4DImageFrom3DImages                     Self;
  typedef itk::ProcessObject                           Superclass;
  typedef itk::SmartPointer< Self >                    Pointer;
  typedef itk::SmartPointer< const Self >              ConstPointer;

  typedef InputImageT                                  InputImageType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( Write4DImageFrom3DImages, ProcessObject );

  /***/
  /***/
  /***/

  /** Set the source image. */
  void SetNumberOfInputImages( unsigned int numInputs );

  /** Set the index'th output image slice */
  void SetNthInputImage( unsigned int outputIndex,
    const InputImageType * img);

  itkSetMacro( FileName, std::string );

  void Write( void )
  { this->Update(); };

protected:
  Write4DImageFrom3DImages( void );
  ~Write4DImageFrom3DImages() {};

  void Update() override;

  void PrintSelf( std::ostream & os, itk::Indent indent ) const override;

private:
  /** itktubeTubeExtractor parameters **/
  Write4DImageFrom3DImages( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) 
    override {};

  typedef itk::Image< typename InputImageType::PixelType, 4 >  OutputImageType;

  unsigned int                       m_NumberOfInputImages;
  typename OutputImageType::Pointer  m_OutputImage;
  std::string                        m_FileName;

};

} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeWrite4DImageFrom3DImages.hxx"
#endif

#endif // End !defined( __tubeWrite4DImageFrom3DImages_h )
