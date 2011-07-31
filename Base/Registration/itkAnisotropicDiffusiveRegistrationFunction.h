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
#ifndef __itkAnisotropicDiffusiveRegistrationFunction_h
#define __itkAnisotropicDiffusiveRegistrationFunction_h

#include "itkPDEDeformableRegistrationFunction.h"

#include "itkAnisotropicDiffusionTensorFunction.h"
#include "itkMeanSquareRegistrationFunction.h"
#include <vector>

namespace itk
{

/** \class itkAnisotropicDiffusiveRegistrationFunction
 * \brief Implements the update function for the
 * itkAnisotropicDiffusiveRegistrationFilter
 *
 * Registration function for registrations using anisotropic diffusive
 * regularizers.  Registration update terms are of the form:
 * Cost = intensityUpdate + div(T1*\grad(u1))v1 + div(T2*\grad(u2))v2 + ...
 * where the types are:
 * - T1..TN are diffusion tensors
 * - u1..uN are deformation vectors
 * - v1..vN are deformation vectors
 * One can specify as many div(T*\grad(u))v regularization terms as you'd like,
 * by passing vectors into ComputeUpdate().
 *
 * Uses the MeanSquareRegistrationFunction as the similarity metric between
 * images.  Uses the implementation in itkAnisotropicDiffusionTensorFunction
 * to calculate the regularization update terms.
 *
 * See: D.F. Pace et al., Deformable image registration of sliding organs using
 * anisotropic diffusive regularization, ISBI 2011.
 *
 * This class is templated over the fixed image type, moving image type and the
 * deformation field type.
 *
 * \sa itkDiffusiveRegistrationFilter
 * \sa itkAnisotropicDiffusiveRegistrationFilter
 * \ingroup FiniteDifferenceFunctions
 * \ingroup Functions
 */

template <class TFixedImage, class TMovingImage, class TDeformationField>
class ITK_EXPORT AnisotropicDiffusiveRegistrationFunction
  : public PDEDeformableRegistrationFunction<TFixedImage,
                                             TMovingImage,
                                             TDeformationField>
{
public:
  /** Standard class typedefs. */
  typedef AnisotropicDiffusiveRegistrationFunction              Self;
  typedef PDEDeformableRegistrationFunction<TFixedImage,
                                            TMovingImage,
                                            TDeformationField>  Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, PDEDeformableRegistrationFunction);

  /** Inherit some parameters from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Convenient typedefs from the superclass. */
  typedef typename Superclass::FixedImageType                FixedImageType;
  typedef typename Superclass::FixedImagePointer             FixedImagePointer;
  typedef typename Superclass::MovingImageType               MovingImageType;
  typedef typename Superclass::MovingImagePointer            MovingImagePointer;
  typedef typename Superclass::DeformationFieldType          DeformationFieldType;
  typedef typename Superclass::DeformationFieldTypePointer   DeformationFieldTypePointer;
  typedef typename Superclass::TimeStepType                  TimeStepType;
  typedef typename Superclass::NeighborhoodType              NeighborhoodType;
  typedef typename Superclass::PixelType                     PixelType;
  typedef typename Superclass::FloatOffsetType               FloatOffsetType;

  /** Deformation field types */
  typedef typename DeformationFieldType::PixelType
      DeformationVectorType;
  typedef typename DeformationVectorType::ValueType
      DeformationVectorComponentType;
  typedef itk::Image< DeformationVectorComponentType, ImageDimension >
      DeformationVectorComponentImageType;
  typedef ZeroFluxNeumannBoundaryCondition< DeformationVectorComponentImageType >
      DeformationVectorComponentImageBoundaryConditionType;
  typedef ConstNeighborhoodIterator
      < DeformationVectorComponentImageType,
      DeformationVectorComponentImageBoundaryConditionType >
      DeformationVectorComponentNeighborhoodType;
  typedef itk::FixedArray
      < DeformationVectorComponentNeighborhoodType, ImageDimension >
       DeformationVectorComponentNeighborhoodArrayType;
  typedef std::vector< DeformationVectorComponentNeighborhoodArrayType >
      DeformationVectorComponentNeighborhoodArrayVectorType;

  /** Typedefs for the intensity-based distance function */
  typedef itk::MeanSquareRegistrationFunction
      < FixedImageType, MovingImageType, DeformationFieldType >
      IntensityDistanceFunctionType;
  typedef typename IntensityDistanceFunctionType::Pointer
      IntensityDistanceFunctionPointer;

  /** Typedefs for the regularization function */
  typedef itk::AnisotropicDiffusionTensorFunction
      < DeformationVectorComponentImageType >
      RegularizationFunctionType;
  typedef typename RegularizationFunctionType::Pointer
      RegularizationFunctionPointer;
  typedef typename RegularizationFunctionType::SpacingType    SpacingType;

  /** Typedefs for the diffusion tensor image */
  typedef typename RegularizationFunctionType::DiffusionTensorType
      DiffusionTensorType;
  typedef typename RegularizationFunctionType::DiffusionTensorImageType
      DiffusionTensorImageType;
  typedef typename RegularizationFunctionType::DiffusionTensorNeighborhoodType
      DiffusionTensorNeighborhoodType;
  typedef std::vector< DiffusionTensorNeighborhoodType >
      DiffusionTensorNeighborhoodVectorType;

  /** Typedefs for the scalar derivatives */
  typedef typename RegularizationFunctionType::ScalarDerivativeType
      ScalarDerivativeType;
  typedef typename RegularizationFunctionType::ScalarDerivativeImageType
      ScalarDerivativeImageType;
  typedef typename RegularizationFunctionType::ScalarDerivativeImageRegionType
      ScalarDerivativeImageRegionType;
  typedef itk::FixedArray < ScalarDerivativeImageRegionType, ImageDimension >
      ScalarDerivativeImageRegionArrayType;
  typedef std::vector < ScalarDerivativeImageRegionArrayType >
      ScalarDerivativeImageRegionArrayVectorType;

  /** Typedefs for the tensor derivatives */
  typedef typename RegularizationFunctionType::TensorDerivativeType
      TensorDerivativeType;
  typedef typename RegularizationFunctionType::TensorDerivativeImageType
      TensorDerivativeImageType;
  typedef typename RegularizationFunctionType::TensorDerivativeImageRegionType
      TensorDerivativeImageRegionType;
  typedef std::vector< TensorDerivativeImageRegionType >
      TensorDerivativeImageRegionVectorType;
  typedef itk::FixedArray < TensorDerivativeImageRegionType, ImageDimension >
      TensorDerivativeImageRegionArrayType;
  typedef std::vector< TensorDerivativeImageRegionArrayType >
      TensorDerivativeImageRegionArrayVectorType;

  /** Typedefs for the multiplication vectors */
  typedef ImageRegionIterator< DeformationFieldType >
      DeformationVectorImageRegionType;
  typedef itk::FixedArray < DeformationVectorImageRegionType, ImageDimension >
      DeformationVectorImageRegionArrayType;
  typedef std::vector< DeformationVectorImageRegionArrayType >
      DeformationVectorImageRegionArrayVectorType;

  /** Computes the time step for an update given a global data structure.
   *  Returns the time step supplied by the user. We don't need
   *  to use the global data supplied since we are returning a fixed value. */
  TimeStepType ComputeGlobalTimeStep(void *itkNotUsed(GlobalData)) const
    { return this->GetTimeStep(); }

  /** Set/Get the time step. For this class of anisotropic diffusion filters,
      the time-step is supplied by the user and remains fixed for all
      updates. */
  void SetTimeStep(const TimeStepType &t)
    {
    m_TimeStep = t;
    if( m_ComputeRegularizationTerm )
      {
      m_RegularizationFunction->SetTimeStep(t);
      }
    // Intensity distance function doesn't have a SetTimeStep(), but it's ok
    // because we only use ComputeUpdate() for it.
    }
  const TimeStepType &GetTimeStep() const
    { return m_TimeStep; }

  /** Utility function to check whether the timestep is stable, optionally based
    * on the spacing of the given image */
  template< class TPixel, unsigned int VImageDimension >
  void CheckTimeStepStability(
      const itk::Image< TPixel, VImageDimension > * input,
      bool useImageSpacing )
    { m_RegularizationFunction->CheckTimeStepStability(input,
                                                       useImageSpacing ); }

  /** Set/get whether to compute the motion field regularization term
   *  Default: true */
  void SetComputeRegularizationTerm( bool compute )
    { m_ComputeRegularizationTerm = compute; }
  bool GetComputeRegularizationTerm() const
    { return m_ComputeRegularizationTerm; }

  /** Set/get whether to compute the intensity distance term
   *  Default: true */
  void SetComputeIntensityDistanceTerm( bool compute )
    { m_ComputeIntensityDistanceTerm = compute; }
  bool GetComputeIntensityDistanceTerm() const
    { return m_ComputeIntensityDistanceTerm; }

  /** Set/get the weighting for the intensity distance update term.  Default
   *  1.0 */
  void SetIntensityDistanceWeighting( double weighting )
    { m_IntensityDistanceWeighting = weighting; }
  double GetIntensityDistanceWeighting() const
    { return m_IntensityDistanceWeighting; }

  /** Set/get the weighting for the regularization update term.  Default 1.0 */
  void SetRegularizationWeighting( double weighting )
    { m_RegularizationWeighting = weighting; }
  double GetRegularizationWeighting() const
    { return m_RegularizationWeighting; }

  /** Set the object's state before each iteration. */
  virtual void InitializeIteration();

  /** Inherited from superclass - do not call this function!  Call the other
   *  ComputeUpdate instead */
  PixelType ComputeUpdate(const NeighborhoodType &neighborhood,
                          void *globalData,
                          const FloatOffsetType &offset = FloatOffsetType(0.0));

  /** Compute the update value. */
  virtual PixelType ComputeUpdate(
      const NeighborhoodType & neighborhood,
      const DiffusionTensorNeighborhoodVectorType & tensorNeighborhoods,
      const ScalarDerivativeImageRegionArrayVectorType
          & intensityFirstDerivatives,
      const TensorDerivativeImageRegionArrayVectorType
          & intensitySecondDerivatives,
      const TensorDerivativeImageRegionVectorType & tensorDerivativeRegions,
      const DeformationVectorImageRegionArrayVectorType
          & multiplicationVectorRegionArrays,
      bool includeInMetricComputations,
      void * globalData,
      const FloatOffsetType& = FloatOffsetType(0.0) );

  /** Returns a pointer to a global data structure that is passed to this
   * object from the solver at each calculation. */
  virtual void * GetGlobalDataPointer() const;

  /** Release the global data structure. */
  virtual void ReleaseGlobalDataPointer(void *GlobalData) const;

  /** Returns the pointers to the regularization function and the intensity
    difference function */
  const RegularizationFunctionType * GetRegularizationFunctionPointer() const
    { return m_RegularizationFunction.GetPointer(); }
  const IntensityDistanceFunctionType * GetIntensityDistanceFunctionPointer()
      const
    { return m_IntensityDistanceFunction.GetPointer(); }

  /** Get the RMS and mean changes in the deformation field. */
  double GetRMSTotalChange() const
    { return this->GetTimeStep() * m_RMSTotalChange; }
  double GetRMSIntensityDistanceChange() const
    { return this->GetTimeStep() * m_RMSIntensityDistanceChange; }
  double GetRMSRegularizationChange() const
    { return this->GetTimeStep() * m_RMSRegularizationChange; };
  double GetMeanTotalChange() const
    { return this->GetTimeStep() * m_MeanTotalChange; }
  double GetMeanIntensityDistanceChange() const
    { return this->GetTimeStep() * m_MeanIntensityDistanceChange; }
  double GetMeanRegularizationChange() const
    { return this->GetTimeStep() * m_MeanRegularizationChange; }

  /** Get the intensity distance and regularization energies, weighted by
   *  the specified weightings. */
  double GetWeightedIntensityDistanceEnergy() const
    { return m_IntensityDistanceWeighting
        * this->GetIntensityDistanceFunctionPointer()->GetEnergy(); }
  double GetWeightedRegularizationEnergy() const
    { return m_RegularizationWeighting * m_RegularizationEnergy; }

protected:
  AnisotropicDiffusiveRegistrationFunction();
  virtual ~AnisotropicDiffusiveRegistrationFunction() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** A global data type for this class of equations.  Used to store information
    for computing the metric and other intermediate products, such as
    derivatives, that may be used by virtual functions called from
    ComputeUpdate().  Caching these values here allows ComputeUpdate() to be
    const and thread-safe. */
  struct GlobalDataStruct
    {
    void *                              m_RegularizationGlobalDataStruct;
    void *                              m_IntensityDistanceGlobalDataStruct;
    double                              m_SumOfSquaredTotalChange;
    double                              m_SumOfSquaredIntensityDistanceChange;
    double                              m_SumOfSquaredRegularizationChange;
    double                              m_SumOfTotalChange;
    double                              m_SumOfIntensityDistanceChange;
    double                              m_SumOfRegularizationChange;
    unsigned long                       m_NumberOfPixelsProcessed;
    };

private:
  // Purposely not implemented
  AnisotropicDiffusiveRegistrationFunction(const Self&);
  void operator=(const Self&);
  // Purposely not implemented

  /** The global timestep. */
  TimeStepType                          m_TimeStep;

  /** The component functions used to calculate the results of this function. */
  RegularizationFunctionPointer         m_RegularizationFunction;
  IntensityDistanceFunctionPointer      m_IntensityDistanceFunction;

  /** Whether or not to compute the intensity distance and motion field
   * regularization terms */
  bool                                  m_ComputeRegularizationTerm;
  bool                                  m_ComputeIntensityDistanceTerm;

  /** Relative weighting between the intensity distance and regularization
   *  terms (default 1,1) */
  double                                m_IntensityDistanceWeighting;
  double                                m_RegularizationWeighting;

  /** Used to calculate the RMS change metrics */
  mutable unsigned long                 m_NumberOfPixelsProcessed;
  mutable double                        m_SumOfSquaredTotalChange;
  mutable double                        m_SumOfSquaredIntensityDistanceChange;
  mutable double                        m_SumOfSquaredRegularizationChange;
  mutable double                        m_RMSTotalChange;
  mutable double                        m_RMSIntensityDistanceChange;
  mutable double                        m_RMSRegularizationChange;
  mutable double                        m_SumOfTotalChange;
  mutable double                        m_SumOfIntensityDistanceChange;
  mutable double                        m_SumOfRegularizationChange;
  mutable double                        m_MeanTotalChange;
  mutable double                        m_MeanIntensityDistanceChange;
  mutable double                        m_MeanRegularizationChange;

  /** Used to calculate the regularization energy. */
  mutable double                        m_RegularizationEnergy;

  /** Mutex locks to protect modifications to metric values. */
  mutable itk::SimpleFastMutexLock      m_MetricCalculationLock;
  mutable itk::SimpleFastMutexLock      m_EnergyCalculationLock;
};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkAnisotropicDiffusiveRegistrationFunction+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkAnisotropicDiffusiveRegistrationFunction.txx"
#endif

#endif
