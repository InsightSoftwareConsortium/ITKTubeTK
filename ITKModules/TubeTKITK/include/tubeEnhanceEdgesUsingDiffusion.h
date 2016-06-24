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
#ifndef __tubeEnhanceEdgesUsingDiffusion_h
#define __tubeEnhanceEdgesUsingDiffusion_h

// ITK includes
#include <itkObject.h>
#include <itktubeAnisotropicEdgeEnhancementDiffusionImageFilter.h>

// TubeTK includes
#include "tubeWrappingMacros.h"

namespace tube
{
/** \class EnhanceEdgesUsingDiffusion
 *
 *  \ingroup TubeTKITK
 */

template< class TInputImage, class TOutputImage >
class EnhanceEdgesUsingDiffusion:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef EnhanceEdgesUsingDiffusion                Self;
  typedef itk::SmartPointer< Self >                 Pointer;
  typedef itk::SmartPointer< const Self >           ConstPointer;

  typedef TInputImage                               InputImageType;
  typedef TOutputImage                              OutputImageType;

  typedef itk::tube::AnisotropicEdgeEnhancementDiffusionImageFilter<
    InputImageType, OutputImageType>  FilterType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( EnhanceEdgesUsingDiffusion, Object );

  /** Set/Get input image */
  tubeWrapSetConstObjectMacro( Input, InputImageType, Filter );
  tubeWrapGetConstObjectMacro( Input, InputImageType, Filter );

  /** Set/Get sigma value */
  tubeWrapSetMacro( Sigma, double, Filter );
  tubeWrapGetMacro( Sigma, double, Filter );

  /** Set/Get ContrastParameterLambdaE value */
  tubeWrapSetMacro( ContrastParameterLambdaE, double, Filter );
  tubeWrapGetMacro( ContrastParameterLambdaE, double, Filter );

  /** Set/Get TimeStep value */
  tubeWrapSetMacro( TimeStep, double, Filter );
  tubeWrapGetMacro( TimeStep, double, Filter );

  /** Set/Get NumberOfIterations value */
  tubeWrapSetMacro( NumberOfIterations, unsigned long, Filter );
  tubeWrapGetMacro( NumberOfIterations, unsigned long, Filter );

  /** Get SigmaOuter value */
  tubeWrapGetMacro( SigmaOuter, double, Filter );

  /** Get ThresholdParameterC value */
  tubeWrapGetMacro( ThresholdParameterC, double, Filter );

  /** Runs the enhancement algorithm */
  tubeWrapUpdateMacro( Filter );

  /** Get output segmentation mask */
  tubeWrapGetObjectMacro( Output, OutputImageType, Filter );

protected:
  EnhanceEdgesUsingDiffusion( void );
  ~EnhanceEdgesUsingDiffusion() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itktubeAnisotropicEdgeEnhancementDiffusionImageFilter parameters **/
  EnhanceEdgesUsingDiffusion( const Self & );
  void operator=( const Self & );

  typename FilterType::Pointer m_Filter;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeEnhanceEdgesUsingDiffusion.hxx"
#endif

#endif // End !defined( __tubeEnhanceEdgesUsingDiffusion_h )
