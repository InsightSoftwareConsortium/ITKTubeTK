/*=========================================================================

Library:   TubeTK/VTree

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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
#ifndef __itkTubeTubeExtractor_h
#define __itkTubeTubeExtractor_h

#include "itkTubeRidgeExtractor.h"
#include "itkTubeRadiusExtractor.h"

namespace itk
{

namespace tube
{

/**
 * This class extract the a tube given an image
 *
 * /sa itkTubeTubeExtractor
 */

template <class TInputImage>
class ITK_EXPORT TubeExtractor : public Object
{
public:

  /**
   * Standard self typedef */
  typedef TubeExtractor             Self;
  typedef Object                    Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /**
   * Run-time type information (and related methods). */
  itkTypeMacro(TubeExtractor, Object);

  itkNewMacro(TubeExtractor);

  /**
   * Type definition for the input image. */
  typedef TInputImage                        ImageType;

  /**
   * Pointer type for the image */
  typedef typename ImageType::Pointer      ImagePointer;

  /**
   * Const Pointer type for the image */
  typedef typename ImageType::ConstPointer ImageConstPointer;

  /**
   * Type definition for the input image pixel type. */
  typedef typename ImageType::PixelType    PixelType;

  /**  Type definition for TubeSpatialObjec*/
  typedef VesselTubeSpatialObject<3>         TubeType;
  typedef typename TubeType::Pointer         TubePointer;
  typedef typename TubeType::TubePointType   TubePointType;

  /**
   * Type definition for the input image pixel type. */
  typedef ContinuousIndex<double, 3>  ContinuousIndexType;

  /**
   * Defines the type of vectors used */
  typedef itk::Vector<double, 3>                     VectorType;

  /**
   * Standard for the number of dimension
   */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      ::itk::GetImageDimension<TInputImage>::ImageDimension);

  /**
   * Set the input image */
  void SetInputImage(ImagePointer inputImage);

  /**
   * Get the input image */
  itkGetConstObjectMacro(Image,ImageType);

  /**
   * Set Data Minimum */
  void SetDataMin(double dataMin);

  /**
   * Get Data Minimum */
  double GetDataMin(void);

  /**
   * Set Data Maximum */
  void SetDataMax(double dataMax);

  /**
   * Get Data Maximum */
  double GetDataMax(void);

  /**
   * Set the radius */
  void SetRadius(double radius);

  /**
   * Get the radius */
  double GetRadius(void);

  /**
   * Set Extract Ridge */
  void ExtractRidge(bool extractRidge);

  /**
   * Get Extract Ridge */
  bool ExtractRidge(void);

  /**
   * Set Extract Valley */
  void ExtractValley(bool extractValley);

  /**
   * Get Extract Valley */
  bool ExtractValley(void);

  /**
   * Get the ridge extractor */
  typename RidgeExtractor<ImageType>::Pointer GetRidgeOp(void);

  /**
   * Get the radius extractor */
  typename RadiusExtractor<ImageType>::Pointer GetRadiusOp(void);

  /**
   * Extract the 3D tube given the position of the first point
   * and the tube ID */
  bool ExtractTube(float x, float y, float z, unsigned int tubeID);

  /**
   * Delete a tube */
  bool DeleteTube(TubeType* tube);

  /**
   * Get the last tube extracted */
  TubePointer GetLastTube(void);

  /**
   * Get the last position */
  ContinuousIndexType   GetLastPosition(void);

  /**
   * Set the tube color */
  void SetColor(float* color);

  /**
   * Get the tube color */
  itkGetMacro(Color,float*);

  /**
   * Set the idle callback */
  void   IdleCallBack(bool (*idleCallBack)());

  /**
   * Set the status callback */
  void   StatusCallBack(void (*statusCallBack)(char *, char *, int));

  /**
   * Set the tube callback */
  void   NewTubeCallBack(void (*newTubeCallBack)(TubeType *));

  /**
   * Set the status callback */
  void   AbortProcess(bool (*abortProcess)());

protected:
  TubeExtractor();
  virtual ~TubeExtractor();
  TubeExtractor(const Self&) {}
  void operator=(const Self&) {}

  typename RidgeExtractor<ImageType>::Pointer  m_RidgeOp;
  typename RadiusExtractor<ImageType>::Pointer m_RadiusOp;
  TubePointer                                  m_Tube;
  bool (*m_IdleCallBack)();
  void (*m_StatusCallBack)(char *, char *, int);
  void (*m_NewTubeCallBack)(TubeType *);
  bool (*m_AbortProcess)();

private:

  bool                     m_Debug;
  ImagePointer             m_Image;
  ContinuousIndexType      m_X;
  float *                  m_Color;

};

} // end namespace tube

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeTubeExtractor.txx"
#endif

#endif /* __itkTubeTubeExtractor_h */
