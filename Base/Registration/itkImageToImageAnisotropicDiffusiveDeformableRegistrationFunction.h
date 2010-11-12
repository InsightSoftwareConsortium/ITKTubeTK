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
#ifndef __itkImageToImageAnisotropicDiffusiveDeformableRegistrationFunction_h
#define __itkImageToImageAnisotropicDiffusiveDeformableRegistrationFunction_h

#include "itkPDEDeformableRegistrationFunction.h"
#include "itkAnisotropicDiffusionTensorFunction.h"
#include "itkMeanSquareRegistrationFunction.h"

namespace itk
{

/** \class itkImageToImageAnisotropicDiffusiveDeformableRegistrationFunction
 * \brief Implements the update function for the
 * itkImageToImageAnisotropicDiffusiveDeformableRegistrationFilter
 *
 * Insert more description + warnings here!!!!
 *
 * This class is templated over the fixed image type, moving image type and the
 * deformation field type.
 *
 * \sa itkImageToImageAnisotropicDiffusiveDeformableRegistrationFilter
 * \ingroup FiniteDifferenceFunctions
 * \ingroup Functions
 */

template <class TFixedImage, class TMovingImage, class TDeformationField>
class ITK_EXPORT ImageToImageAnisotropicDiffusiveDeformableRegistrationFunction
  : public PDEDeformableRegistrationFunction<TFixedImage,
                                             TMovingImage,
                                             TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef ImageToImageAnisotropicDiffusiveDeformableRegistrationFunction
      Self;
  typedef PDEDeformableRegistrationFunction<TFixedImage,
                                            TMovingImage,
                                            TDeformationField>  Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, PDEDeformableRegistrationFunction);

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Types for the moving image. */
  typedef typename Superclass::MovingImageType        MovingImageType;
  typedef typename Superclass::MovingImagePointer     MovingImagePointer;

  /** Types for the fixed image. */
  typedef typename Superclass::FixedImageType         FixedImageType;
  typedef typename Superclass::FixedImagePointer      FixedImagePointer;

  /** Types for the deformation field. */
  typedef typename Superclass::DeformationFieldType   DeformationFieldType;
  typedef typename Superclass::DeformationFieldTypePointer
      DeformationFieldTypePointer;

  /** Inherit some types from the superclass. */
  typedef typename Superclass::TimeStepType           TimeStepType;
  typedef typename Superclass::NeighborhoodType       NeighborhoodType;
  typedef typename Superclass::PixelType              PixelType;
  typedef typename Superclass::FloatOffsetType        FloatOffsetType;

  /** Deformation field types - types for the deformation vectors, deformation
   *  vector components, and vector component images
   */
  typedef typename DeformationFieldType::PixelType    DeformationVectorType;
  typedef typename DeformationVectorType::ValueType
      DeformationVectorComponentType;
  typedef itk::Image< DeformationVectorComponentType, ImageDimension >
      DeformationVectorComponentImageType;

  /** Normal vector types */
  typedef double                                   NormalVectorComponentType;
  typedef itk::Vector< NormalVectorComponentType, ImageDimension >
      NormalVectorType;
  typedef itk::Image< NormalVectorType, ImageDimension >
      NormalVectorImageType;

  /** Typedefs for the intensity-based distance function */
  typedef itk::MeanSquareRegistrationFunction< FixedImageType,
                                               MovingImageType,
                                               DeformationFieldType >
                                               IntensityDistanceFunctionType;
  typedef typename IntensityDistanceFunctionType::Pointer
      IntensityDistanceFunctionPointer;

  /** Typedefs for the regularization function */
  typedef itk::AnisotropicDiffusionTensorFunction
      < DeformationVectorComponentImageType >
      RegularizationFunctionType;
  typedef typename RegularizationFunctionType::Pointer
      RegularizationFunctionPointer;

  /** Typedefs for the diffusion tensor image */
  typedef typename RegularizationFunctionType::DiffusionTensorImageType
      DiffusionTensorImageType;
  typedef typename DiffusionTensorImageType::Pointer
      DiffusionTensorImagePointer;

  /** Boundary condition typedefs (defined in RegularizationFunction for
    * diffusion tensors) */
  typedef ZeroFluxNeumannBoundaryCondition< DeformationVectorComponentImageType >
      DeformationVectorComponentImageBoundaryConditionType;
  typedef ZeroFluxNeumannBoundaryCondition< NormalVectorImageType >
      NormalVectorImageBoundaryConditionType;

  /** Neighborhood iterator typedefs */
  typedef ConstNeighborhoodIterator
      < DeformationVectorComponentImageType,
      DeformationVectorComponentImageBoundaryConditionType >
      DeformationVectorComponentNeighborhoodIteratorType;
  typedef itk::FixedArray
      < DeformationVectorComponentNeighborhoodIteratorType, ImageDimension >
       DeformationVectorComponentNeighborhoodIteratorArrayType;
  typedef ConstNeighborhoodIterator
      < NormalVectorImageType, NormalVectorImageBoundaryConditionType >
      NormalVectorImageNeighborhoodIteratorType;

  typedef typename RegularizationFunctionType::DiffusionTensorNeighborhoodType
      DiffusionTensorNeighborhoodIteratorType;

  /** Set/Get the time step for an update */
  void SetTimeStep(const TimeStepType &t)
    {
    m_TimeStep = t;
    m_RegularizationFunction->SetTimeStep(t);
    // Intensity distance function doesn't have a SetTimeStep(), but it's ok
    // because we only use ComputeUpdate() for that function
    }
  const TimeStepType &GetTimeStep() const
    { return m_TimeStep; }
  /** Computes the time step for an update given a global data structure.  For
   * this class of anisotropic diffusion filters, the time-step is supplied
   * by the user and is fixed for all updates, so the global data structure
   * is not used.
   */
  virtual TimeStepType ComputeGlobalTimeStep(void * itkNotUsed(GlobalData) )
      const
    { return m_TimeStep; }

  /** Set/get whether to compute the motion field regularization term
   *  Default: true
   */
  void SetComputeRegularizationTerm( bool compute )
    { m_ComputeRegularizationTerm = compute; }
  bool GetComputeRegularizationTerm() const
    { return m_ComputeRegularizationTerm; }

  /** Set/get whether to compute the intensity distance term
   *  Default: true
   */
  void SetComputeIntensityDistanceTerm( bool compute )
    { m_ComputeIntensityDistanceTerm = compute; }
  bool GetComputeIntensityDistanceTerm() const
    { return m_ComputeIntensityDistanceTerm; }

  /** Set the object's state before each iteration. */
  virtual void InitializeIteration();

  /** Inherited from superclass - do not call this function!  Call the other
   *  ComputeUpdate instead
   */
  PixelType ComputeUpdate(const NeighborhoodType & neighborhood,
                          void * globalData,
                          const FloatOffsetType & offset );

  /** Compute the update value. */
  virtual PixelType ComputeUpdate(
                     const NeighborhoodType &neighborhood,
                     const NormalVectorImageNeighborhoodIteratorType
                              &normalVectorImageNeighborhood,
                     const DiffusionTensorNeighborhoodIteratorType
                              &tangentialNeighborhoodTensor,
                     const DeformationVectorComponentNeighborhoodIteratorArrayType
                              &tangentialNeighborhoodDeformationFieldComponents,
                     const DiffusionTensorNeighborhoodIteratorType
                              &normalNeighborhoodTensor,
                     const DeformationVectorComponentNeighborhoodIteratorArrayType
                              &normalNeighborhoodDeformationFieldComponents,
                     void *globalData,
                     const FloatOffsetType& = FloatOffsetType(0.0));

  /** Returns a pointer to a global data structure that is passed to this
   * object from the solver at each calculation.*/
  virtual void *GetGlobalDataPointer() const;

  /** Release the global data structure. */
  virtual void ReleaseGlobalDataPointer(void *GlobalData) const;

protected:
  ImageToImageAnisotropicDiffusiveDeformableRegistrationFunction();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** A global data type for this class of equations.  Used to store information
    for computing the metric and other intermediate products, such as
    derivatives, that may be used by virtual functions called from
    ComputeUpdate().  Caching these values here allows ComputeUpdate() to be
    const and thread-safe.*/
  struct GlobalDataStruct
    {
      void *                            m_RegularizationGlobalDataStruct;
      void *                            m_IntensityDistanceGlobalDataStruct;
    };

private:
  // Purposely not implemented
  ImageToImageAnisotropicDiffusiveDeformableRegistrationFunction(const Self&);
  void operator=(const Self&); // Purposely not implemented

  /** The global timestep. */
  TimeStepType                          m_TimeStep;

  /** The component functions used to calculate the results of this function. */
  RegularizationFunctionPointer         m_RegularizationFunction;
  IntensityDistanceFunctionPointer      m_IntensityDistanceFunction;

  /** Mutex lock to protect modification to metric. */
  mutable SimpleFastMutexLock           m_MetricCalculationLock;

  /** Whether or not to compute the intensity distance and motion field
   * regularization terms */
  bool                                  m_ComputeRegularizationTerm;
  bool                                  m_ComputeIntensityDistanceTerm;

};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkImageToImageAnisotropicDiffusiveDeformableRegistrationFunction+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkImageToImageAnisotropicDiffusiveDeformableRegistrationFunction.txx"
#endif

#endif
