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

#ifndef itktubeMeanSquareRegistrationFunction_h
#define itktubeMeanSquareRegistrationFunction_h

#include <itkCentralDifferenceImageFunction.h>
#include <itkCovariantVector.h>
#include <itkInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkPDEDeformableRegistrationFunction.h>
#include <itkPoint.h>
#include <mutex>

namespace itk
{

namespace tube
{

/** \class MeanSquareRegistrationFunction
 *
 * This class encapsulate the PDE which drives the demons registration
 * algorithm. It is used by MeanSquareRegistrationFilter to compute the
 * output deformation field which will map a moving image onto a
 * a fixed image.
 *
 * Non-integer moving image values are obtained by using
 * interpolation. The default interpolator is of type
 * LinearInterpolateImageFunction. The user may set other
 * interpolators via method SetMovingImageInterpolator. Note that the input
 * interpolator must derive from base class InterpolateImageFunction.
 *
 * This class is templated over the fixed image type, moving image type,
 * and the deformation field type.
 *
 * \warning This filter assumes that the fixed image type, moving image type
 * and deformation field type all have the same number of dimensions.
 *
 * \sa MeanSquareRegistrationFilter
 * \ingroup FiniteDifferenceFunctions
 */
template <class TFixedImage, class TMovingImage, class TDeformationField>
class MeanSquareRegistrationFunction
  : public PDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
{
public:
  /** Standard class type alias. */
  using Self = MeanSquareRegistrationFunction;
  using Superclass = PDEDeformableRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro(MeanSquareRegistrationFunction);

  /** MovingImage image type. */
  using MovingImageType = typename Superclass::MovingImageType;
  using MovingImagePointer = typename Superclass::MovingImagePointer;
  using MovingImagePixelType = typename MovingImageType::PixelType;

  /** FixedImage image type. */
  using FixedImageType = typename Superclass::FixedImageType;
  using FixedImagePointer = typename Superclass::FixedImagePointer;
  using IndexType = typename FixedImageType::IndexType;
  using SizeType = typename FixedImageType::SizeType;
  using SpacingType = typename FixedImageType::SpacingType;

  /** Deformation field type. */
  typedef typename Superclass::DisplacementFieldType DeformationFieldType;
  typedef typename DeformationFieldType::Pointer     DeformationFieldPointer;
  typedef typename DeformationFieldType::PixelType   DeformationFieldPixelType;

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Inherit some enums from the superclass. */
  using PixelType = typename Superclass::PixelType;
  using RadiusType = typename Superclass::RadiusType;
  using NeighborhoodType = typename Superclass::NeighborhoodType;
  using FloatOffsetType = typename Superclass::FloatOffsetType;
  using TimeStepType = typename Superclass::TimeStepType;

  /** Interpolator type. */
  using CoordinateType = double;
  using InterpolatorType = InterpolateImageFunction<MovingImageType, CoordinateType>;
  using InterpolatorPointer = typename InterpolatorType::Pointer;
  using PointType = typename InterpolatorType::PointType;
  using DefaultInterpolatorType = LinearInterpolateImageFunction<MovingImageType, CoordinateType>;

  /** Covariant vector type. */
  using CovariantVectorType = CovariantVector<double, itkGetStaticConstMacro(ImageDimension)>;

  /** Gradient calculator type. */
  using GradientCalculatorType = CentralDifferenceImageFunction<FixedImageType>;
  typedef typename GradientCalculatorType::Pointer GradientCalculatorPointer;

  /** Set the moving image interpolator. */
  void
  SetMovingImageInterpolator(InterpolatorType * ptr)
  {
    m_MovingImageInterpolator = ptr;
  }

  /** Get the moving image interpolator. */
  InterpolatorType *
  GetMovingImageInterpolator(void)
  {
    return m_MovingImageInterpolator;
  }

  /** This class uses a constant time step of 1. */
  virtual TimeStepType
  ComputeGlobalTimeStep(void * itkNotUsed(globalData)) const override
  {
    return m_TimeStep;
  }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  virtual void *
  GetGlobalDataPointer(void) const override
  {
    GlobalDataStruct * global = new GlobalDataStruct();
    return global;
  }

  /** Release memory for global data structure. */
  virtual void
  ReleaseGlobalDataPointer(void * GlobalData) const override
  {
    delete (GlobalDataStruct *)GlobalData;
  }

  /** Set the object's state before each iteration. */
  virtual void
  InitializeIteration(void) override;

  /** This method is called by a finite difference solver image filter at
   * each pixel that does not lie on a data set boundary */
  virtual PixelType
  ComputeUpdate(const NeighborhoodType & neighborhood,
                void *                   globalData,
                const FloatOffsetType &  offset = FloatOffsetType(0.0)) override;

  /** Computes the intensity difference between the fixed and moving image
   *  at the given index, under the given deformation vector. */
  virtual double
  ComputeIntensityDifference(const IndexType & index, const DeformationFieldPixelType & itvec);

  /** Get the energy mutex lock  */
  void
  SetEnergy(double energy)
  {
#ifndef wasi__
    std::lock_guard<std::mutex> mutexHolder(m_EnergyCalculationLock);
#endif
    this->m_Energy = energy;
  }

  void
  SetBackgroundIntensity(MovingImagePixelType intensity)
  {
    m_BackgroundIntensity = intensity;
  }
  MovingImagePixelType
  GetBackgroundIntensity(void) const
  {
    return m_BackgroundIntensity;
  }

  void
  SetIntensityDifferenceThreshold(double threshold)
  {
    m_IntensityDifferenceThreshold = threshold;
  }
  double
  GetIntensityDifferenceThreshold(void) const
  {
    return m_IntensityDifferenceThreshold;
  }


protected:
  MeanSquareRegistrationFunction(void);
  ~MeanSquareRegistrationFunction(void) {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  /** FixedImage image neighborhood iterator type. */
  using FixedImageNeighborhoodIteratorType = ConstNeighborhoodIterator<FixedImageType>;

  /** A global data type for this class of equation. Used to store
   * iterators for the fixed image. */
  struct GlobalDataStruct
  {
    FixedImageNeighborhoodIteratorType m_FixedImageIterator;
  }; // End struct GlobalDataStruct

private:
  MeanSquareRegistrationFunction(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  /** Cache fixed image information. */
  SpacingType m_FixedImageSpacing;

  /** Function to compute derivatives of the fixed image. */
  GradientCalculatorPointer m_FixedImageGradientCalculator;

  /** Function to interpolate the moving image. */
  InterpolatorPointer m_MovingImageInterpolator;

  /** The global time step. */
  TimeStepType m_TimeStep;

  /** Threshold below which the denominator term is considered zero. */
  double m_DenominatorThreshold;

  /** Threshold below which two intensity value are assumed to match. */
  double m_IntensityDifferenceThreshold;

#ifndef wasi__
  mutable std::mutex m_EnergyCalculationLock;
#endif

  MovingImagePixelType m_BackgroundIntensity;

}; // End class MeanSquareRegistrationFunction

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeMeanSquareRegistrationFunction.hxx"
#endif

#endif // End !defined( __itktubeMeanSquareRegistrationFunction_h )
