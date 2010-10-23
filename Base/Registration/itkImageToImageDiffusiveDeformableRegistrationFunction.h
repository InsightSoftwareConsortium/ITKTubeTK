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
#ifndef __itkImageToImageDiffusiveDeformableRegistrationFunction_h
#define __itkImageToImageDiffusiveDeformableRegistrationFunction_h

#include "itkPDEDeformableRegistrationFunction.h"
#include "itkAnisotropicDiffusionTensorFunction.h"
#include "itkMeanSquareRegistrationFunction.h"

namespace itk
{

/** \class itkImageToImageDiffusiveDeformableRegistrationFunction
 * \brief Implements the update function for the
 * itkImageToImageDiffusiveDeformableRegistrationFilter
 *
 * Insert more description + warnings here!!!!
 *
 * This class is templated over the fixed image type, moving image type and the
 * deformation field type.
 *
 * \sa itkImageToImageDiffusiveDeformableRegistrationFilter
 * \ingroup FiniteDifferenceFunctions
 * \ingroup Functions
 */

template <class TFixedImage, class TMovingImage, class TDeformationField>
class ITK_EXPORT ImageToImageDiffusiveDeformableRegistrationFunction
  : public PDEDeformableRegistrationFunction<TFixedImage,
                                             TMovingImage,
                                             TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef ImageToImageDiffusiveDeformableRegistrationFunction   Self;
  typedef PDEDeformableRegistrationFunction<TFixedImage,
                                            TMovingImage,
                                            TDeformationField>  Superclass;
  typedef SmartPointer<Self>                                    Pointer;
  typedef SmartPointer<const Self>                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, PDEDeformableRegistrationFunction);

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType        MovingImageType;
  typedef typename Superclass::MovingImagePointer     MovingImagePointer;

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType         FixedImageType;
  typedef typename Superclass::FixedImagePointer      FixedImagePointer;
  typedef typename FixedImageType::IndexType          IndexType;
  typedef typename FixedImageType::SizeType           SizeType;
  typedef typename FixedImageType::SpacingType        SpacingType;

  /** Deformation field type. */
  typedef typename Superclass::DeformationFieldType
                                            DeformationFieldType;
  typedef typename Superclass::DeformationFieldTypePointer
                                            DeformationFieldTypePointer;
  typedef typename Superclass::DeformationFieldType::ConstPointer
                                            DeformationFieldConstPointer;

  /** Inherit some enums from the superclass. */
  typedef typename Superclass::TimeStepType           TimeStepType;
  typedef typename Superclass::RadiusType             RadiusType;
  typedef typename Superclass::NeighborhoodType       NeighborhoodType;
  typedef typename Superclass::PixelType              PixelType;
  typedef typename Superclass::FloatOffsetType        FloatOffsetType;

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Typedefs for the deformation field components */
  // ex. vector < double, 3 >
  typedef typename DeformationFieldType::PixelType
                                DeformationFieldVectorType;
  // ex. double
  typedef typename DeformationFieldType::PixelType::ValueType
                                DeformationFieldScalarType;
  // ex. image of doubles
  typedef itk::Image< DeformationFieldScalarType, ImageDimension >
                                DeformationFieldComponentImageType;
  typedef typename DeformationFieldComponentImageType::ConstPointer
                                DeformationFieldComponentImageConstPointer;

  typedef ZeroFluxNeumannBoundaryCondition< DeformationFieldComponentImageType >
                                DeformationFieldComponentBoundaryConditionType;
  typedef ConstNeighborhoodIterator<DeformationFieldComponentImageType,
                                DeformationFieldComponentBoundaryConditionType>
                                DeformationFieldComponentNeighborhoodType;
  typedef itk::FixedArray< DeformationFieldComponentNeighborhoodType, ImageDimension >
                                DeformationFieldComponentNeighborhoodArrayType;

  /** Typedefs for normal vector */
  typedef DeformationFieldVectorType                     NormalVectorType;
  typedef itk::Image< NormalVectorType, ImageDimension > NormalVectorImageType;

  typedef ZeroFluxNeumannBoundaryCondition< NormalVectorImageType >
                                NormalVectorImageBoundaryConditionType;
  typedef ConstNeighborhoodIterator< NormalVectorImageType,
                                NormalVectorImageBoundaryConditionType >
                                NormalVectorImageNeighborhoodType;

  /** Typedefs for the regularization function used by this function */
  typedef itk::AnisotropicDiffusionTensorFunction
                                          < DeformationFieldComponentImageType >
                                          RegularizationFunctionType;
  typedef typename RegularizationFunctionType::Pointer
                                          RegularizationFunctionPointer;

  /** Typedefs for the diffusion tensor image */
  typedef typename RegularizationFunctionType::DiffusionTensorImageType
                                          DiffusionTensorImageType;
  typedef typename DiffusionTensorImageType::Pointer
                                          DiffusionTensorImagePointer;
  typedef typename RegularizationFunctionType::DefaultBoundaryConditionType
                                          DefaultBoundaryConditionType;
  typedef typename RegularizationFunctionType::DiffusionTensorNeighborhoodType
                                          DiffusionTensorNeighborhoodType;

  /** Typedefs for the FiniteDifferenceFunction intensity-based distance
    * function used by this function */
  typedef itk::MeanSquareRegistrationFunction< FixedImageType,
                                               MovingImageType,
                                               DeformationFieldType >
                                               IntensityDistanceFunctionType;
  typedef typename IntensityDistanceFunctionType::Pointer
                                          IntensityDistanceFunctionPointer;

  /** Set the object's state before each iteration. */
  virtual void InitializeIteration();

  /** Inherited from superclass - do not call this function!  Call ComputeUpdate
   * (neighborhood, normalVectorImageNeighborhood, tangentialNeighborhoodTensor,
   * tangentialNeighborhoodDeformationFieldComponents, normalNeighborhoodTensor,
   * normalNeighborhoodDeformationFieldComponents, globalData, offset) instead
   */
  virtual PixelType ComputeUpdate(const NeighborhoodType &neighborhood,
                                  void *globalData,
                                  const FloatOffsetType &offset
                                                        = FloatOffsetType(0.0));

  /** Compute the equation value. */
  virtual PixelType ComputeUpdate(
                     const NeighborhoodType &neighborhood,
                     const NormalVectorImageNeighborhoodType
                              &normalVectorImageNeighborhood,
                     const DiffusionTensorNeighborhoodType
                              &tangentialNeighborhoodTensor,
                     const DeformationFieldComponentNeighborhoodArrayType
                              &tangentialNeighborhoodDeformationFieldComponents,
                     const DiffusionTensorNeighborhoodType
                              &normalNeighborhoodTensor,
                     const DeformationFieldComponentNeighborhoodArrayType
                              &normalNeighborhoodDeformationFieldComponents,
                     void *globalData,
                     const FloatOffsetType& = FloatOffsetType(0.0));

  /** This class uses a constant timestep of 1. */
  virtual TimeStepType ComputeGlobalTimeStep(void * itkNotUsed( GlobalData ))
      const
    {
    return m_TimeStep;
    }

  /** Set/Get the time step. For this class of anisotropic diffusion filters,
      the time-step is supplied by the user and remains fixed for all
      updates. */
  void SetTimeStep(const TimeStepType &t)
    {
    m_TimeStep = t;
    m_RegularizationFunction->SetTimeStep(t);
    // Intensity distance function doesn't have a SetTimeStep(), but it's ok
    // because we only use ComputeUpdate() from that function
    }

  /** This class uses a constant timestep of 1. */
  const TimeStepType &GetTimeStep() const
    {
    return m_TimeStep;
    }

  /** Whether to compute the motion field regularization term (for testing)
   *  Default: true */
  void SetComputeRegularizationTerm( bool compute );
  bool GetComputeRegularizationTerm() const;

  /** Whether to compute the intensity distance term (for testing)
   *  Default: true */
  void SetComputeIntensityDistanceTerm( bool compute );
  bool GetComputeIntensityDistanceTerm() const;

  /** Returns a pointer to a global data structure that is passed to this
   * object from the solver at each calculation.*/
  virtual void *GetGlobalDataPointer() const;

  /** Release the global data structure. */
  virtual void ReleaseGlobalDataPointer(void *GlobalData) const;

protected:
  ImageToImageDiffusiveDeformableRegistrationFunction();
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
  ImageToImageDiffusiveDeformableRegistrationFunction(const Self&);
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
# include "Templates/itkImageToImageDiffusiveDeformableRegistrationFunction+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkImageToImageDiffusiveDeformableRegistrationFunction.txx"
#endif

#endif
