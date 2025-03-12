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

#ifndef itktubeAnisotropicDiffusiveRegistrationFunction_h
#define itktubeAnisotropicDiffusiveRegistrationFunction_h

#include "itktubeAnisotropicDiffusionTensorFunction.h"
#include "itktubeMeanSquareRegistrationFunction.h"

#include <itkPDEDeformableRegistrationFunction.h>

#include <vector>
#ifndef wasi__
#  include <mutex>
#endif

namespace itk
{

namespace tube
{

/** \class AnisotropicDiffusiveRegistrationFunction
 * \brief Implements the update function for the
 * itktubeAnisotropicDiffusiveRegistrationFilter
 *
 * Registration function for registrations using anisotropic diffusive
 * regularizers.  Registration update terms are of the form:
 * Cost = intensityUpdate + div( T1*\grad( u1 ) )v1
 * + div( T2*\grad( u2 ) )v2 + ...
 * where the types are:
 * - T1..TN are diffusion tensors
 * - u1..uN are deformation vectors
 * - v1..vN are deformation vectors
 * One can specify as many div( T*\grad( u ) )v regularization terms as
 * you'd like,
 * by passing vectors into ComputeUpdate().
 *
 * Uses the MeanSquareRegistrationFunction as the similarity metric between
 * images.  Uses the implementation in AnisotropicDiffusionTensorFunction
 * to calculate the regularization update terms.
 *
 * See: D.F. Pace et al., Deformable image registration of sliding organs
 * using
 * anisotropic diffusive regularization, ISBI 2011.
 *
 * This class is templated over the fixed image type, moving image type and
 * the
 * deformation field type.
 *
 * \sa DiffusiveRegistrationFilter
 * \sa AnisotropicDiffusiveRegistrationFilter
 * \ingroup FiniteDifferenceFunctions
 * \ingroup Functions
 */

template <class TFixedImage, class TMovingImage, class TDeformationField>
class AnisotropicDiffusiveRegistrationFunction
  : public PDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
{
public:
  /** Standard class type alias. */
  using Self = AnisotropicDiffusiveRegistrationFunction;
  using Superclass = PDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ). */
  itkTypeMacro(AnisotropicDiffusiveRegistrationFunction, PDEDeformableRegistrationFunction);

  /** Inherit some parameters from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Convenient type alias from the superclass. */
  using FixedImageType = typename Superclass::FixedImageType;
  using FixedImagePointer = typename Superclass::FixedImagePointer;
  using MovingImageType = typename Superclass::MovingImageType;
  using MovingImagePointer = typename Superclass::MovingImagePointer;
  using MovingImagePixelType = typename MovingImageType::PixelType;
  using DeformationFieldType = typename Superclass::DisplacementFieldType;
  using DeformationFieldPointer = typename DeformationFieldType::Pointer;
  using TimeStepType = typename Superclass::TimeStepType;
  using NeighborhoodType = typename Superclass::NeighborhoodType;
  using PixelType = typename Superclass::PixelType;
  using FloatOffsetType = typename Superclass::FloatOffsetType;

  /** Deformation field types */
  typedef typename DeformationFieldType::PixelType  DeformationVectorType;
  typedef typename DeformationVectorType::ValueType DeformationVectorComponentType;
  using DeformationVectorComponentImageType = itk::Image<DeformationVectorComponentType, ImageDimension>;
  using DeformationVectorComponentImageBoundaryConditionType =
    ZeroFluxNeumannBoundaryCondition<DeformationVectorComponentImageType>;
  using DeformationVectorComponentNeighborhoodType =
    ConstNeighborhoodIterator<DeformationVectorComponentImageType,
                              DeformationVectorComponentImageBoundaryConditionType>;
  using DeformationVectorComponentNeighborhoodArrayType =
    itk::FixedArray<DeformationVectorComponentNeighborhoodType, ImageDimension>;
  using DeformationVectorComponentNeighborhoodArrayVectorType =
    std::vector<DeformationVectorComponentNeighborhoodArrayType>;

  /** Typedefs for the intensity-based distance function */
  using IntensityDistanceFunctionType =
    MeanSquareRegistrationFunction<FixedImageType, MovingImageType, DeformationFieldType>;
  typedef typename IntensityDistanceFunctionType::Pointer IntensityDistanceFunctionPointer;

  /** Typedefs for the regularization function */
  using RegularizationFunctionType = AnisotropicDiffusionTensorFunction<DeformationVectorComponentImageType>;
  typedef typename RegularizationFunctionType::Pointer     RegularizationFunctionPointer;
  typedef typename RegularizationFunctionType::SpacingType SpacingType;

  /** Typedefs for the diffusion tensor image */
  typedef typename RegularizationFunctionType::DiffusionTensorType             DiffusionTensorType;
  typedef typename RegularizationFunctionType::DiffusionTensorImageType        DiffusionTensorImageType;
  typedef typename RegularizationFunctionType::DiffusionTensorNeighborhoodType DiffusionTensorNeighborhoodType;
  using DiffusionTensorNeighborhoodVectorType = std::vector<DiffusionTensorNeighborhoodType>;

  /** Typedefs for the scalar derivatives */
  typedef typename RegularizationFunctionType::ScalarDerivativeType            ScalarDerivativeType;
  typedef typename RegularizationFunctionType::ScalarDerivativeImageType       ScalarDerivativeImageType;
  typedef typename RegularizationFunctionType::ScalarDerivativeImageRegionType ScalarDerivativeImageRegionType;
  using ScalarDerivativeImageRegionArrayType = itk::FixedArray<ScalarDerivativeImageRegionType, ImageDimension>;
  using ScalarDerivativeImageRegionArrayVectorType = std::vector<ScalarDerivativeImageRegionArrayType>;

  /** Typedefs for the tensor derivatives */
  typedef typename RegularizationFunctionType::TensorDerivativeType            TensorDerivativeType;
  typedef typename RegularizationFunctionType::TensorDerivativeImageType       TensorDerivativeImageType;
  typedef typename RegularizationFunctionType::TensorDerivativeImageRegionType TensorDerivativeImageRegionType;
  using TensorDerivativeImageRegionVectorType = std::vector<TensorDerivativeImageRegionType>;
  using TensorDerivativeImageRegionArrayType = itk::FixedArray<TensorDerivativeImageRegionType, ImageDimension>;
  using TensorDerivativeImageRegionArrayVectorType = std::vector<TensorDerivativeImageRegionArrayType>;

  /** Typedefs for the multiplication vectors */
  using DeformationVectorImageRegionType = ImageRegionIterator<DeformationFieldType>;
  using DeformationVectorImageRegionArrayType = itk::FixedArray<DeformationVectorImageRegionType, ImageDimension>;
  using DeformationVectorImageRegionArrayVectorType = std::vector<DeformationVectorImageRegionArrayType>;

  /** Computes the time step for an update given a global data structure.
   *  Returns the time step supplied by the user. We don't need
   *  to use the global data supplied since we are returning a fixed
   *  value. */
  TimeStepType
  ComputeGlobalTimeStep(void * itkNotUsed(globalData)) const override
  {
    return this->GetTimeStep();
  }

  /** Set/Get the time step. For this class of anisotropic diffusion
   * filters, the time-step is supplied by the user and remains fixed for
   * all updates. */
  void
  SetTimeStep(const TimeStepType & t)
  {
    m_TimeStep = t;
    if (m_ComputeRegularizationTerm)
    {
      m_RegularizationFunction->SetTimeStep(t);
    }
    // Intensity distance function doesn't have a SetTimeStep(), but it's ok
    // because we only use ComputeUpdate() for it.
  }
  const TimeStepType &
  GetTimeStep(void) const
  {
    return m_TimeStep;
  }

  /** Utility function to check whether the timestep is stable, optionally
   * based on the spacing of the given image */
  template <class TPixel, unsigned int VImageDimension>
  void
  CheckTimeStepStability(const itk::Image<TPixel, VImageDimension> * input, bool useImageSpacing)
  {
    m_RegularizationFunction->CheckTimeStepStability(input, useImageSpacing);
  }

  /** Set/get whether to compute the motion field regularization term
   *  Default: true */
  void
  SetComputeRegularizationTerm(bool compute)
  {
    m_ComputeRegularizationTerm = compute;
  }
  bool
  GetComputeRegularizationTerm(void) const
  {
    return m_ComputeRegularizationTerm;
  }

  /** Set/get whether to compute the intensity distance term
   *  Default: true */
  void
  SetComputeIntensityDistanceTerm(bool compute)
  {
    m_ComputeIntensityDistanceTerm = compute;
  }
  bool
  GetComputeIntensityDistanceTerm(void) const
  {
    return m_ComputeIntensityDistanceTerm;
  }

  /** Set/get the weighting for the regularization update term.  Default
   * 1.0 */
  void
  SetRegularizationWeighting(double weighting)
  {
    m_RegularizationWeighting = weighting;
  }
  double
  GetRegularizationWeighting(void) const
  {
    return m_RegularizationWeighting;
  }

  /** Set/get the background intensity in the moving image.  Default 0.0 */
  void
  SetBackgroundIntensity(MovingImagePixelType bg)
  {
    m_IntensityDistanceFunction->SetBackgroundIntensity(bg);
  }

  MovingImagePixelType
  GetBackgroundIntensity(void) const
  {
    return m_IntensityDistanceFunction->GetBackgroundIntensity();
  }

  /** Set the object's state before each iteration. */
  virtual void
  InitializeIteration(void) override;

  /** Inherited from superclass - do not call this function!  Call the other
   *  ComputeUpdate instead */
  PixelType
  ComputeUpdate(const NeighborhoodType & neighborhood,
                void *                   globalData,
                const FloatOffsetType &  offset = FloatOffsetType(0.0)) override;

  /** Compute the update value.  The intensityDistanceTerm and
   *  regularizationTerm are outputs.  Incorporates weighting between
   *  intensity distance term and regularization term, but does not yet
   *  incorporate the time step. */
  virtual PixelType
  ComputeUpdate(const NeighborhoodType &                            neighborhood,
                const DiffusionTensorNeighborhoodVectorType &       tensorNeighborhoods,
                const ScalarDerivativeImageRegionArrayVectorType &  deformationComponentFirstOrderDerivativeRegions,
                const TensorDerivativeImageRegionArrayVectorType &  deformationComponentSecondOrderDerivativeRegions,
                const TensorDerivativeImageRegionVectorType &       tensorDerivativeRegions,
                const DeformationVectorImageRegionArrayVectorType & multiplicationVectorRegionArrays,
                void *                                              globalData,
                PixelType &                                         intensityDistanceTerm,
                PixelType &                                         regularizationTerm,
                const FloatOffsetType & = FloatOffsetType(0.0));

  /** Updates the energy associated with the intensity distance term */
  virtual double
  ComputeIntensityDistanceEnergy(const typename NeighborhoodType::IndexType index,
                                 const DeformationVectorType &              update);

  /** Updates the energy associated with the regularization */
  virtual double
  ComputeRegularizationEnergy(
    const DiffusionTensorNeighborhoodVectorType &      tensorNeighborhoods,
    const ScalarDerivativeImageRegionArrayVectorType & deformationComponentFirstOrderDerivativeRegions);

  /** Returns a pointer to a global data structure that is passed to this
   * object from the solver at each calculation. */
  virtual void *
  GetGlobalDataPointer(void) const override;

  /** Release the global data structure. */
  virtual void
  ReleaseGlobalDataPointer(void * GlobalData) const override;

  /** Returns the pointers to the regularization function and the intensity
    difference function */
  const RegularizationFunctionType *
  GetRegularizationFunctionPointer(void) const
  {
    return m_RegularizationFunction.GetPointer();
  }
  const IntensityDistanceFunctionType *
  GetIntensityDistanceFunctionPointer(void) const
  {
    return m_IntensityDistanceFunction.GetPointer();
  }

  /** Get the intensity distance and regularization energies.
   *  The regularization energy incorporates the weighting between the
   *  intensity distance term and the regularization term. */
  double
  GetIntensityDistanceEnergy(void) const
  {
    return this->GetIntensityDistanceFunctionPointer()->GetEnergy();
  }
  double
  GetRegularizationEnergy(void) const
  {
    return m_RegularizationEnergy;
  }

protected:
  AnisotropicDiffusiveRegistrationFunction(void);
  virtual ~AnisotropicDiffusiveRegistrationFunction(void) {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** Returns the update from the regularization component */
  virtual PixelType
  ComputeRegularizationUpdate(
    const DiffusionTensorNeighborhoodVectorType &       tensorNeighborhoods,
    const ScalarDerivativeImageRegionArrayVectorType &  deformationComponentFirstOrderDerivativeRegions,
    const TensorDerivativeImageRegionArrayVectorType &  deformationComponentSecondOrderDerivativeRegions,
    const TensorDerivativeImageRegionVectorType &       tensorDerivativeRegions,
    const DeformationVectorImageRegionArrayVectorType & multiplicationVectorRegionArrays,
    void *                                              globalData,
    const FloatOffsetType & = FloatOffsetType(0.0));

  //  /** Update the RMS and mean update change statistics */
  //  virtual void UpdateRMSAndMeanUpdateStatistics(
  //    const PixelType & updateTerm,
  //    const PixelType & intensityDistanceTerm,
  //    const PixelType & regularizationTerm,
  //    const TimeStepType & timestep,
  //    void * globalData );

  /** A global data type for this class of equations.  Used to store
   * information for computing the metric and other intermediate products,
   * such as derivatives, that may be used by virtual functions called from
   * ComputeUpdate().  Caching these values here allows ComputeUpdate() to
   * be const and thread-safe. */
  struct GlobalDataStruct
  {
    void * m_RegularizationGlobalDataStruct;
    void * m_IntensityDistanceGlobalDataStruct;
  }; // End struct GlobalDataStruct

private:
  // Purposely not implemented
  AnisotropicDiffusiveRegistrationFunction(const Self &);
  void
  operator=(const Self &);
  // Purposely not implemented

  /** The global time step. */
  TimeStepType m_TimeStep;

  /** The component functions used to calculate the results of this
   * function. */
  RegularizationFunctionPointer    m_RegularizationFunction;
  IntensityDistanceFunctionPointer m_IntensityDistanceFunction;

  /** Whether or not to compute the intensity distance and motion field
   * regularization terms */
  bool m_ComputeRegularizationTerm;
  bool m_ComputeIntensityDistanceTerm;

  /** Relative weighting between the intensity distance and regularization
   *  terms ( default 1 ) */
  double m_RegularizationWeighting;

  /** Used to calculate the regularization energy. */
  mutable double m_RegularizationEnergy;

  /** Mutex locks to protect modifications to metric values. */
#ifndef wasi__
  mutable std::mutex m_MetricCalculationLock;
  mutable std::mutex m_EnergyCalculationLock;
#endif

}; // End class AnisotropicDiffusiveRegistrationFunction

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeAnisotropicDiffusiveRegistrationFunction.hxx"
#endif

// End !defined( __itktubeAnisotropicDiffusiveRegistrationFunction_h )
#endif
