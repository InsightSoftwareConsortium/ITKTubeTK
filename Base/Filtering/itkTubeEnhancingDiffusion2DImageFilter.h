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
#ifndef __itkTubeEnhancingDiffusion2DImageFilter_h
#define __itkTubeEnhancingDiffusion2DImageFilter_h

#include "itkImageToImageFilter.h"
#include <vector>

namespace itk
{

/** \class TubeEnhancingDiffusion2DImageFilter
 *
 * Complete rewrite of previous versions, only using itk/vnl routines
 * for derivatives, eigensystem calculations and diffusion. Internally,
 * the input image image is converted to internalprecision (float) for
 * calculation, and converted back when returning the results.
 *
 * Uses simple forward Euler scheme (explicit) with 3x3 stencil,
 * see eg phd of Joachim Weickert for theory and implementation regarding
 * the construction of this discretization scheme. See 'Tube Enhancing
 * Diffusion', Manniesing, media 2006, for information regarding the
 * construction of the diffusion tensor.
 *
 * - Stores all elements of the Hessian of the complete image during
 *   diffusion. An alternative implementation is to only store the
 *   scale for which the vesselness has maximum response, and to
 *   recalculate the hessian (locally) during diffusion. Also stores
 *   the current image, ie at iteration i + temp image, therefore the complete
 *   memory consumption approximately peaks at 8 times the input image
 *   (input image in float)
 * - The hessian is stored as six individual images, an alternative
 *   implementation is to use the itk symmetric second rank tensor
 *   as pixeltype (and eg using the class SymmetricEigenAnalysisImage
 *   Filter). However, we are lazy, and using this since we rely
 *   on vnl datatypes and its eigensystem calculations
 * - note: most of computation time is spent at calculation of vesselness
 *   response
 *
 * - PixelT         short, 2D
 *   Precision      float, 2D
 *
 *
 * - todo
 *   - using parallelism/threading eg over scales
 *   - completely itk-fying, eg eigenvalues calculation
 *   - possibly embedding within itk-diffusion framework
 *   - itk expert to have a look at use of iterators
 *     (there must be a potential gain there)
 *
 * email: r.manniesing@erasmusmc.nl
 *
 */
template <class PixelT = short int, unsigned int DimensionT = 2>
class ITK_EXPORT TubeEnhancingDiffusion2DImageFilter :
    public ImageToImageFilter<Image<PixelT, DimensionT> ,
                              Image<PixelT, DimensionT> >
{

public:

  typedef float                                           Precision;
  typedef Image<PixelT, DimensionT>                       ImageType;
  typedef Image<Precision, DimensionT>                    PrecisionImageType;

  typedef TubeEnhancingDiffusion2DImageFilter             Self;
  typedef ImageToImageFilter<ImageType,ImageType>         Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  itkNewMacro(Self);
  itkTypeMacro(TubeEnhancingDiffusion2DImageFilter, ImageToImageFilter);

  itkSetMacro(TimeStep, Precision);
  itkSetMacro(Iterations, unsigned int);
  itkSetMacro(RecalculateTubeness, unsigned int);

  itkSetMacro(Beta, Precision);
  itkSetMacro(Gamma, Precision);

  itkSetMacro(Epsilon, Precision);
  itkSetMacro(Omega, Precision);
  itkSetMacro(Sensitivity, Precision);

  void SetScales(const std::vector<Precision> scales)
    {
    m_Scales = scales;
    }

  itkBooleanMacro(DarkObjectLightBackground);
  itkSetMacro(DarkObjectLightBackground,bool);
  itkBooleanMacro(Verbose);
  itkSetMacro(Verbose,bool);

  // some defaults for lowdose example
  // used in the paper
  void SetDefaultPars()
    {
    m_TimeStep                  = 0.25;
    m_Iterations                = 200;
    m_RecalculateTubeness     = 100;
    m_Beta                      = 0.5;
    m_Gamma                     = 5.0;
    m_Epsilon                   = 0.01;
    m_Omega                     = 25.0;
    m_Sensitivity               = 20.0;

    m_Scales.resize(2);
    m_Scales[0] = 6;
    m_Scales[1] = 8;

    m_DarkObjectLightBackground = true;
    m_Verbose                   = true;
    }

protected:
  TubeEnhancingDiffusion2DImageFilter();
  ~TubeEnhancingDiffusion2DImageFilter() {};
  void PrintSelf(std::ostream &os, Indent indent) const;
  void GenerateData();

private:

  TubeEnhancingDiffusion2DImageFilter(const Self&);
  void operator=(const Self&);

  Precision                 m_TimeStep;
  unsigned int              m_Iterations;
  unsigned int              m_RecalculateTubeness;
  Precision                 m_Beta;
  Precision                 m_Gamma;
  Precision                 m_Epsilon;
  Precision                 m_Omega;
  Precision                 m_Sensitivity;
  std::vector<Precision>    m_Scales;
  bool                      m_DarkObjectLightBackground;
  bool                      m_Verbose;

  unsigned int              m_CurrentIteration;

  // current hessian for which we have max vesselresponse
  typename PrecisionImageType::Pointer m_Dxx;
  typename PrecisionImageType::Pointer m_Dxy;
  typename PrecisionImageType::Pointer m_Dyy;

  void VED2DSingleIteration( typename PrecisionImageType::Pointer );

  // Calculates maxvessel response of the range
  // of scales and stores the hessian of each voxel
  // into the member images m_Dij.
  void MaxTubeResponse( const typename PrecisionImageType::Pointer );

  // calculates diffusion tensor
  // based on current values of hessian (for which we have
  // maximim vessel response).
  void DiffusionTensor();

  // Sorted increasing magnitude: l1, l2
  inline Precision TubenessFunction2D ( const Precision, const Precision );

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeEnhancingDiffusion2DImageFilter.txx"
#endif

#endif
