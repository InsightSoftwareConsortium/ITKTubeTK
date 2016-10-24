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
#ifndef __tubeEnhanceTubesUsingDiffusion_h
#define __tubeEnhanceTubesUsingDiffusion_h

// ITK includes
#include <itkProcessObject.h>

// TubeTK includes
#include "tubeWrappingMacros.h"

#include "itktubeTubeEnhancingDiffusion2DImageFilter.h"

namespace tube
{
/** \class EnhanceTubesUsingDiffusion
 *
 *  \ingroup TubeTKITK
 *
 *  Compute vesselness score of an image using Frangi's method.
 */

template< class TPixel, unsigned int Dimension >
class EnhanceTubesUsingDiffusion:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef EnhanceTubesUsingDiffusion                Self;
  typedef itk::ProcessObject                        Superclass;
  typedef itk::SmartPointer< Self >                 Pointer;
  typedef itk::SmartPointer< const Self >           ConstPointer;

  typedef itk::Image< TPixel, Dimension >           ImageType;

  typedef itk::tube::TubeEnhancingDiffusion2DImageFilter< TPixel,
    Dimension >                                     FilterType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(EnhanceTubesUsingDiffusion, ProcessObject);

  /** Set/Get minimum sigma/scale */
  itkSetMacro(MinSigma, double);
  itkGetConstMacro(MinSigma, double);

  /** Set/Get maximum sigma/scale */
  itkSetMacro(MaxSigma, double);
  itkGetConstMacro(MaxSigma, double);

  /** Set/Get number of sigma/scale steps */
  itkSetMacro(NumSigmaSteps, unsigned int);
  itkGetConstMacro(NumSigmaSteps, unsigned int);

  /** Set/Get for how many iterations do we recalculate tubeness */
  tubeWrapSetMacro(RecalculateTubeness, unsigned int, Filter);
  tubeWrapGetMacro(RecalculateTubeness, unsigned int, Filter);

  /** Set/Get sensitive of the filter to blobness. default = 0.5 */
  tubeWrapSetMacro(Beta, double, Filter);
  tubeWrapGetMacro(Beta, double, Filter);

  /** Set/Get sensitive of the filter to second order structureness.
   *  default = 5.0 */
  tubeWrapSetMacro(Gamma, double, Filter);
  tubeWrapGetMacro(Gamma, double, Filter);

  /** Set/Get epsilon. default = 0.01 */
  tubeWrapSetMacro(Epsilon, double, Filter);
  tubeWrapGetMacro(Epsilon, double, Filter);

  /** Set/Get Omega. default = 25.0 */
  tubeWrapSetMacro(Omega, double, Filter);
  tubeWrapGetMacro(Omega, double, Filter);

  /** Set/Get Sensitivity. default = 20.0 */
  tubeWrapSetMacro(Sensitivity, double, Filter);
  tubeWrapGetMacro(Sensitivity, double, Filter);

  /** Set/Get time step */
  tubeWrapSetMacro(TimeStep, double, Filter);
  tubeWrapGetMacro(TimeStep, double, Filter);

  /** Set/Get iterations */
  tubeWrapSetMacro(Iterations, unsigned int, Filter);
  tubeWrapGetMacro(Iterations, unsigned int, Filter);

  /** Set/Get input image */
  tubeWrapSetConstObjectMacro(Input, ImageType, Filter);
  tubeWrapGetConstObjectMacro(Input, ImageType, Filter);

  /** Compute vesselness image */
  void Update( void );

  /** Get output vesselness image */
  tubeWrapGetObjectMacro(Output, ImageType, Filter);

protected:
  EnhanceTubesUsingDiffusion( void );
  ~EnhanceTubesUsingDiffusion() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itkEnhanceTubesUsingDiffusionFilter parameters **/
  EnhanceTubesUsingDiffusion(const Self &);
  void operator=(const Self &);

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) {};

  double                          m_MinSigma;
  double                          m_MaxSigma;
  unsigned int                    m_NumSigmaSteps;
  typename FilterType::Pointer    m_Filter;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeEnhanceTubesUsingDiffusion.hxx"
#endif

#endif // End !defined( __tubeEnhanceTubesUsingDiffusion_h )
