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
#ifndef __itkTubeTubeNetworkExtractor_h
#define __itkTubeTubeNetworkExtractor_h

#include "itkTubeTubeExtractor.h"

namespace itk
{

namespace tube
{

/**
 * This class extract the a tube given an image
 *
 * /sa itkTubeTubeNetworkExtractor
 */

template <class TInputImage, class TInputMask>
class ITK_EXPORT TubeNetworkExtractor : public TubeExtractor<TInputImage>
{
public:

  /**
   * Standard self typedef */
  typedef TubeNetworkExtractor             Self;
  typedef Object                           Superclass;
  typedef SmartPointer<Self>               Pointer;
  typedef SmartPointer<const Self>         ConstPointer;

  itkTypeMacro( Self, Object );

  itkNewMacro( Self );

  /**
   * Type definition for the input image. */
  typedef TInputImage                      ImageType;

  /**
   * Type definition for the input mask. */
  typedef TInputMask                       MaskType;

  /**
   * Input image pixel type. */
  typedef typename TInputImage::PixelType  PixelType;

  /**
   * Input image index type. */
  typedef typename TInputImage::IndexType  IndexType;

  /**
   * Standard for the number of dimension
   */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Typedef for VesselTubeSpatialObject */
  typedef VesselTubeSpatialObject< ImageDimension >    TubeType;

  /**
   * Defines the type of vectors used
   */
  typedef Vector<double, ImageDimension>               VectorType;

  /**
   * Set the input image */
  void SetInputImage( typename ImageType::Pointer inputImage );

  /**
   * Get the input image */
  itkGetConstObjectMacro( Image, ImageType );

  /**
   * Get the tube Network */
  typename TubeType::Pointer GetTubeNetwork( void );

  /**
   * Set use mask */
  itkSetMacro( AEUseMask,bool );

  /**
   * Get use mask */
  itkGetMacro( AEUseMask,bool );

  /**
   * Set the mask  */
  void SetAutoExtractMask( typename MaskType::Pointer autoExtractMask );

  /**
   * Get the mask  */
  itkGetConstObjectMacro( AEMask, MaskType );

  /**
   * Auto extract tubes using a mask */
  double  AutoExtractThresh( void );

  /**
   * Auto extract tubes using a mask */
  void   AutoExtractThresh( double newAEThresh );

  /**
   * AutoExtract tubes */
  bool   AutoExtract( int zMin, int zMax );

  /**
   * Generate output mask */
  void   DrawMask( MaskType * maskImage );

  /**
   * Save a tube net */
  //bool   Save( char *fname );

  /**
   * Load a tube net*/
  //bool   Load( char *fname );

  /**
   * Set the tube callback */
  void   NewTubeCallBack( void ( *newTubeCallBack )( TubeType * ) );

protected:

  TubeNetworkExtractor( void );
  virtual ~TubeNetworkExtractor( void );
  TubeNetworkExtractor( const Self& ) {}
  void operator=( const Self& ) {}

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  typename ImageType::Pointer      m_Image;
  typename TubeType::Pointer       m_TubeNetwork;
  int                              m_TubeNum;

  bool                             m_AEUseMask;
  typename MaskType::Pointer       m_AEMask;
  float                            m_AEThresh;

}; // End class TubeNetworkExtractor

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeTubeNetworkExtractor.txx"
#endif

#endif // End !defined(__itkTubeTubeNetworkExtractor_h)
