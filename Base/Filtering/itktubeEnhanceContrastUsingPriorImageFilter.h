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

#ifndef __itktubeEnhanceContrastUsingPriorImageFilter_h
#define __itktubeEnhanceContrastUsingPriorImageFilter_h

#include <itkArray.h>
#include <itkImageToImageFilter.h>
#include <itktubeContrastCostFunction.h>
#include <itkNormalVariateGenerator.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkFRPROptimizer.h>
#include <itkImageRegionIterator.h>

namespace itk
{

namespace tube
{

/** \class EnhanceContrastUsingPriorImageFilter
 */

template< class TPixel, unsigned int VDimension >
class EnhanceContrastUsingPriorImageFilter
  : public ImageToImageFilter< Image< TPixel, VDimension >,
  Image< TPixel, VDimension > >
{
public:
  /** Standard class typedefs. */
  typedef Image< TPixel, VDimension >                     ImageType;
  typedef EnhanceContrastUsingPriorImageFilter            Self;
  typedef ImageToImageFilter< ImageType,ImageType >       Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( EnhanceContrastUsingPriorImageFilter, ImageToImageFilter );

  /** Some convenient typedefs. */
  typedef itk::tube::ContrastCostFunction< TPixel, VDimension >
                                                ContrastCostFunctionType;
  typedef itk::OnePlusOneEvolutionaryOptimizer  InitialOptimizerType;
  typedef itk::FRPROptimizer                    OptimizerType;
  typedef itk::ImageRegionIterator< ImageType > ImageIteratorType;

  /** Set/Get input Mask Image */
  itkSetObjectMacro( InputMaskImage, ImageType );
  itkGetObjectMacro( InputMaskImage, ImageType );

  /** Set/Get Object Scale */
  itkSetMacro( ObjectScale, float );
  itkGetMacro( ObjectScale, float );

  /** Set/Get Background Scale */
  itkSetMacro( BackgroundScale, float );
  itkGetMacro( BackgroundScale, float );

  /** Set/Get Mask Object Value */
  itkSetMacro( MaskObjectValue, int );
  itkGetMacro( MaskObjectValue, int );

  /** Set/Get Mask Background Value */
  itkSetMacro( MaskBackgroundValue, int );
  itkGetMacro( MaskBackgroundValue, int );

  /** Set/Get Optimization Iterations */
  itkSetMacro( OptimizationIterations, int );
  itkGetMacro( OptimizationIterations, int );

  /** Set/Get Optimization Seed */
  itkSetMacro( OptimizationSeed, int );
  itkGetMacro( OptimizationSeed, int );

protected:
  EnhanceContrastUsingPriorImageFilter( void );
  virtual ~EnhanceContrastUsingPriorImageFilter( void ) {}

  void PrintSelf(std::ostream& os, Indent indent) const;
  virtual void GenerateData( void );

private:
  EnhanceContrastUsingPriorImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  typename ImageType::Pointer            m_InputMaskImage;
  float                                  m_ObjectScale;
  float                                  m_BackgroundScale;
  int                                    m_MaskObjectValue;
  int                                    m_MaskBackgroundValue;
  int                                    m_OptimizationIterations;
  int                                    m_OptimizationSeed;

}; // End class EnhanceContrastUsingPriorImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeEnhanceContrastUsingPriorImageFilter.hxx"
#endif

#endif // End !defined(__itktubeEnhanceContrastUsingPriorImageFilter_h)
