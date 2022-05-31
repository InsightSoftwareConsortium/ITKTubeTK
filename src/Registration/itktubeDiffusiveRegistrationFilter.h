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

#ifndef __itktubeDiffusiveRegistrationFilter_h
#define __itktubeDiffusiveRegistrationFilter_h

#include "itktubeAnisotropicDiffusiveRegistrationFunction.h"

#include <itkPDEDeformableRegistrationFilter.h>
#include <itkVectorIndexSelectionCastImageFilter.h>

namespace itk
{

namespace tube
{

struct EnergiesStruct
{
  double  TotalEnergy;
  double  IntensityDistanceEnergy;
  double  RegularizationEnergy;

  void zero( void )
    {
    TotalEnergy = 0.0;
    IntensityDistanceEnergy = 0.0;
    RegularizationEnergy = 0.0;
    }

  void copyFrom( const EnergiesStruct & rhs )
    {
    TotalEnergy = rhs.TotalEnergy;
    IntensityDistanceEnergy = rhs.IntensityDistanceEnergy;
    RegularizationEnergy = rhs.RegularizationEnergy;
    }

  void difference( const EnergiesStruct & lhs, const EnergiesStruct & rhs )
    {
    TotalEnergy = lhs.TotalEnergy - rhs.TotalEnergy;
    IntensityDistanceEnergy = lhs.IntensityDistanceEnergy
      - rhs.IntensityDistanceEnergy;
    RegularizationEnergy = lhs.RegularizationEnergy
      - rhs.RegularizationEnergy;
    }
}; // End struct EnergiesStruct

/** \class DiffusiveRegistrationFilter
 * \brief Registration filter for registrations using anisotropic diffusive
 * regularizers, for example for sliding organ registration.
 *
 * Traditional deformable image registration imposes a uniform
 * smoothness constraint on the deformation field. This
 * is not appropriate when registering images visualizing organs
 * that slide relative to each other, and therefore leads to registration
 * inaccuracies.
 *
 * This filter is a base class for non-parametric deformable registration
 * algorithms with regularization terms based on anisotropic diffusion.  For
 * example, these regularizers can accomodate deformation field
 * discontinuities
 * that are expected when considering sliding motion.
 *
 * The update term is composed of two parts: an intensity distance term and
 * a
 * regularization term.  The intensity distance term uses the sum of square
 * difference metric, so this registration algorithm is appropriate for
 * monomodal image registration term only.
 *
 * The update term for the regularization will be of the form:
 * div( T1*\grad( u1 ) )v1 + div( T2*\grad( u2 ) )v2 + ...
 * + div( TN*\grad( uN ) )vN
 * where the types are:
 * - T1..TN are diffusion tensors
 * - u1..uN are deformation vectors
 * - v1..vN are deformation vectors
 * It is assumed that T's and v's are constant throughout the registration,
 * and
 * so can be precomputed once, while u's must be updated on each
 * registration iteration.
 *
 * This base class implements the diffusive regularization.  Algorithms
 * implementing anisotropic regularization should derive it and override the
 * following functions:
 * - GetNumberOfTerms(): returns the number of div( T*\grad( u ) )v terms
 * - ComputeDiffusionTensorImages(): allocate and populate the T images
 * - InitializeDeformationComponentAndDerivativeImages(): allocate the u
 *   images
 * and their derivatives
 * - ComputeMultiplicationVectorImages(): allocate and populate the v images
 * - UpdateDeformationComponentImages(): update the u images at each
 *   iteration
 * See itktubeAnisotropicDiffusiveRegistrationFilter for an example derived
 * filter.
 *
 * See: D.F. Pace et al., Deformable image registration of sliding organs
 * using
 * anisotropic diffusive regularization, ISBI 2011.
 *
 * This class is templated over the type of the fixed image, the type of the
 * moving image and the type of the deformation field.
 *
 * \sa AnisotropicDiffusiveRegistrationFunction
 * \sa AnisotropicDiffusiveRegistrationFilter
 * \ingroup DeformableImageRegistration
 * \ingroup MultiThreaded
 */

template< class TFixedImage, class TMovingImage, class TDeformationField >
class DiffusiveRegistrationFilter
  : public PDEDeformableRegistrationFilter< TFixedImage, TMovingImage,
                                            TDeformationField >
{
public:
  /** Standard class typedefs. */
  typedef DiffusiveRegistrationFilter                         Self;
  typedef PDEDeformableRegistrationFilter< TFixedImage,
                                          TMovingImage,
                                          TDeformationField > Superclass;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;

  /** Method for creation through the object factory.  Usually defined with
    * itkNewMacro(), but the type of the registration function depends on
    * the
    * type of this object.  Can't call the overridden function
    * CreateRegistrationFunction() from the base class constructor, so we'll
    * call it here. Derived classes should use this instead of
    * itkNewMacro(). */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( DiffusiveRegistrationFilter,
    PDEDeformableRegistrationFilter );

  /** Inherit some parameters from the superclass. */
  itkStaticConstMacro( ImageDimension, unsigned int,
    Superclass::ImageDimension );

  /** Convenient typedefs from the superclass. */
  typedef typename Superclass::FixedImageType        FixedImageType;
  typedef typename Superclass::FixedImagePointer     FixedImagePointer;
  typedef typename Superclass::MovingImageType       MovingImageType;
  typedef typename Superclass::MovingImagePointer    MovingImagePointer;
  typedef typename MovingImageType::PixelType        MovingImagePixelType;
  typedef typename Superclass::DisplacementFieldType DeformationFieldType;
  typedef typename Superclass::TimeStepType          TimeStepType;

  typedef typename Superclass::DisplacementFieldPointer
      DeformationFieldPointer;
  typedef typename Superclass::FiniteDifferenceFunctionType
      FiniteDifferenceFunctionType;

  /** Typedefs used in multithreading */
  typedef typename Superclass::OutputImageType          OutputImageType;
  typedef typename Superclass::OutputImagePointer       OutputImagePointer;
  typedef typename Superclass::UpdateBufferType         UpdateBufferType;
  typedef typename UpdateBufferType::RegionType         ThreadRegionType;

  /** Output image and update buffer types */
  typedef itk::ImageRegionIterator< OutputImageType >
      OutputImageRegionType;
  typedef typename FiniteDifferenceFunctionType::NeighborhoodType
      NeighborhoodType;

  /** The registration function type */
  typedef AnisotropicDiffusiveRegistrationFunction
      < FixedImageType, MovingImageType, DeformationFieldType >
      RegistrationFunctionType;
  typedef typename RegistrationFunctionType::RegularizationFunctionType
      RegularizationFunctionType;
  typedef typename RegistrationFunctionType::RegularizationFunctionPointer
      RegularizationFunctionPointer;
  typedef typename RegistrationFunctionType::SpacingType SpacingType;

  /** Deformation component types ( i.e. component of a deformation field,
   *  still a vector */
  typedef std::vector< DeformationFieldPointer >
      DeformationFieldArrayType;
  typedef typename RegistrationFunctionType::DeformationVectorType
      DeformationVectorType;

  /** Deformation vector component types ( i.e. scalar within a vector ) */
  typedef typename RegistrationFunctionType::DeformationVectorComponentType
      DeformationVectorComponentType;
  typedef typename RegistrationFunctionType
    ::DeformationVectorComponentImageType
      DeformationVectorComponentImageType;
  typedef typename DeformationVectorComponentImageType::Pointer
      DeformationVectorComponentImagePointer;
  typedef typename itk::FixedArray< DeformationVectorComponentImagePointer,
    ImageDimension >
      DeformationComponentImageArrayType;
  typedef typename RegistrationFunctionType
    ::DeformationVectorComponentNeighborhoodType
      DeformationVectorComponentNeighborhoodType;
  typedef typename DeformationVectorComponentImageType::RegionType
      ThreadDeformationVectorComponentImageRegionType;

  /** Diffusion tensor image types */
  typedef typename RegistrationFunctionType::DiffusionTensorType
      DiffusionTensorType;
  typedef typename RegistrationFunctionType::DiffusionTensorImageType
      DiffusionTensorImageType;
  typedef typename DiffusionTensorImageType::Pointer
      DiffusionTensorImagePointer;
  typedef std::vector< DiffusionTensorImagePointer >
      DiffusionTensorImageArrayType;
  typedef typename RegistrationFunctionType::DiffusionTensorNeighborhoodType
      DiffusionTensorNeighborhoodType;
  typedef typename
      RegistrationFunctionType::DiffusionTensorNeighborhoodVectorType
      DiffusionTensorNeighborhoodVectorType;
  typedef typename DiffusionTensorImageType::RegionType
      ThreadDiffusionTensorImageRegionType;

  /** Scalar derivative image types */
  typedef typename RegistrationFunctionType::ScalarDerivativeImageType
      ScalarDerivativeImageType;
  typedef typename ScalarDerivativeImageType::Pointer
      ScalarDerivativeImagePointer;
  typedef typename itk::FixedArray< ScalarDerivativeImagePointer >
      ScalarDerivativeImageArrayType;
  typedef std::vector< ScalarDerivativeImageArrayType >
      ScalarDerivativeImageArrayVectorType;
  typedef typename RegistrationFunctionType::ScalarDerivativeImageRegionType
      ScalarDerivativeImageRegionType;
  typedef typename
      RegistrationFunctionType::ScalarDerivativeImageRegionArrayVectorType
      ScalarDerivativeImageRegionArrayVectorType;
  typedef typename ScalarDerivativeImageType::RegionType
      ThreadScalarDerivativeImageRegionType;

  /** Tensor derivative matrix image types */
  typedef typename RegistrationFunctionType::TensorDerivativeImageType
      TensorDerivativeImageType;
  typedef typename TensorDerivativeImageType::Pointer
      TensorDerivativeImagePointer;
  typedef std::vector< TensorDerivativeImagePointer >
      TensorDerivativeImageVectorType;
  typedef typename itk::FixedArray< TensorDerivativeImagePointer >
      TensorDerivativeImageArrayType;
  typedef std::vector< TensorDerivativeImageArrayType >
      TensorDerivativeImageArrayVectorType;
  typedef typename RegistrationFunctionType::TensorDerivativeImageRegionType
      TensorDerivativeImageRegionType;
  typedef typename
      RegistrationFunctionType::TensorDerivativeImageRegionVectorType
      TensorDerivativeImageRegionVectorType;
  typedef typename
      RegistrationFunctionType::TensorDerivativeImageRegionArrayVectorType
      TensorDerivativeImageRegionArrayVectorType;
  typedef typename TensorDerivativeImageType::RegionType
      ThreadTensorDerivativeImageRegionType;

  /** Typedefs for the multiplication vectors */
  typedef typename itk::FixedArray< DeformationFieldPointer >
      DeformationVectorImageArrayType;
  typedef std::vector< DeformationVectorImageArrayType >
      DeformationVectorImageArrayVectorType;
  typedef typename RegistrationFunctionType
    ::DeformationVectorImageRegionType
      DeformationVectorImageRegionType;
  typedef typename
      RegistrationFunctionType::DeformationVectorImageRegionArrayVectorType
      DeformationVectorImageRegionArrayVectorType;

  /** Stopping criterion mask types */
  typedef FixedImageType    StoppingCriterionMaskImageType;
  typedef FixedImagePointer StoppingCriterionMaskPointer;
  typedef typename StoppingCriterionMaskImageType::RegionType
    ThreadStoppingCriterionMaskImageRegionType;
  typedef ImageRegionIterator< StoppingCriterionMaskImageType >
    StoppingCriterionMaskImageRegionType;

  /** Convenience functions to set/get the registration functions
   * timestep. */
  virtual void SetTimeStep( const TimeStepType t )
    { m_OriginalTimeStep = t; }
  virtual const TimeStepType& GetTimeStep( void ) const
    { return m_OriginalTimeStep; }

  /** Set/get whether to compute the motion field regularization term
   *  Default: true */
  void SetComputeRegularizationTerm( bool compute )
    { this->GetRegistrationFunctionPointer()->
      SetComputeRegularizationTerm( compute ); }
  bool GetComputeRegularizationTerm( void ) const
    { return this->GetRegistrationFunctionPointer()->
      GetComputeRegularizationTerm(); }

  /** Set/get whether to compute the intensity distance term
   *  Default: true */
  void SetComputeIntensityDistanceTerm( bool compute )
    { this->GetRegistrationFunctionPointer()->
      SetComputeIntensityDistanceTerm( compute ); }
  bool GetComputeIntensityDistanceTerm( void ) const
    { return this->GetRegistrationFunctionPointer()->
      GetComputeIntensityDistanceTerm(); }

  /** Set/get the weightings for the regularization update term.  If
   *  using multiresolution registration and the current level is past the
   *  length of the weight vector, the last weight in the vector will be
   *  used.
   *  Default: 1.0 */
  void SetRegularizationWeightings( std::vector< double >& weightings )
    { m_RegularizationWeightings = weightings; }
  const std::vector< double >& GetRegularizationWeightings( void ) const
    { return m_RegularizationWeightings; }

  /** Set/get the background intensity of the moving image, used by the
   *  intensity distance function.  Default 0.0 */
  void SetBackgroundIntensity( MovingImagePixelType bg )
    {
    this->GetRegistrationFunctionPointer()->SetBackgroundIntensity( bg );
    }
  MovingImagePixelType GetBackgroundIntensity( void ) const
    {
    return this->GetRegistrationFunctionPointer()->GetBackgroundIntensity();
    }

  /** The number of div( T\grad( u ) )v terms we sum for the regularizer.
   *  Reimplement in derived classes. */
  virtual int GetNumberOfTerms( void ) const
    { return 1; }

  /** Set/get a pointer to an image that is to be used for the template when
   *  computing member images.  This is usually the original fixed image.
   *  The
   *  attributes of this filter's output are used if the high resolution
   *  template is not set.  For proper behavior, you must set this if using
   *  a multiresolution registration. */
  virtual void SetHighResolutionTemplate( FixedImageType * templateImage )
    { m_HighResolutionTemplate = templateImage; }
  virtual FixedImageType * GetHighResolutionTemplate( void )
    { return m_HighResolutionTemplate; }

  /** Get current resolution level being processed. */
  itkGetConstReferenceMacro( CurrentLevel, unsigned int );

  /** Set/get a mask in which the RMS error does not contribute to the
   * stopping criterion.  Any non-zero voxels will not be considered when
   * determining the stopping criterion. */
  void SetStoppingCriterionMask( StoppingCriterionMaskImageType * mask )
    { m_StoppingCriterionMask = mask; }
  StoppingCriterionMaskImageType * GetStoppingCriterionMask( void ) const
    { return m_StoppingCriterionMask; }

  /** Set/get the number of iterations that elapse between evaluations of
   * the stopping criterion.  Default 50. */
  void SetStoppingCriterionEvaluationPeriod( unsigned int numIterations )
    { m_StoppingCriterionEvaluationPeriod = numIterations; }
  unsigned int GetStoppingCriterionEvaluationPeriod( void ) const
    { return m_StoppingCriterionEvaluationPeriod; }

  /** Set/get the total energy change, summed over the stopping criterion
   * evaluation
   *  period, that indicates the stopping criterion has been reached.
   *  Default 0,
   *  indicating that the total energy is not used and the stopping
   *  criterion is
   *  defined solely by the number of iterations specified or the maximum
   *  RMS change.
   *  Specify in absolute value. */
  void SetStoppingCriterionMaxTotalEnergyChange( double energyChange )
    { m_StoppingCriterionMaxTotalEnergyChange = energyChange; }
  double GetStoppingCriterionMaxTotalEnergyChange( void ) const
    { return m_StoppingCriterionMaxTotalEnergyChange; }

protected:
  DiffusiveRegistrationFilter( void );
  virtual ~DiffusiveRegistrationFilter( void ) {}

  void PrintSelf( std::ostream& os, Indent indent ) const override;

  /** Handy for array indexing. */
  enum DivTerm { GAUSSIAN };

  /** Initialization occuring before the registration iterations begin. */
  virtual void Initialize( void ) override;

  /** Allocate images used during the registration. */
  virtual void AllocateImageMembers( void );

  /** Allocate the deformation component images and their derivative images.
   *  ( which may be updated throughout the registration ). Reimplement in
   *  derived
   *  classes. */
  virtual void InitializeDeformationComponentAndDerivativeImages( void );

  /** Allocate and populate the diffusion tensor images.
   *  Reimplement in derived classes. */
  virtual void ComputeDiffusionTensorImages( void );

  /** Computes the first-order partial derivatives of the diffusion tensor
   *  images.  Override in derived classes if the diffusion tensor image
   *  pointers are not unique, to avoid computing the derivatives multiple
   *  times. */
  virtual void ComputeDiffusionTensorDerivativeImages( void );

  /** Helper to compute the first-order partial derivatives of the diffusion
   *  tensor images */
  virtual void ComputeDiffusionTensorDerivativeImageHelper(
      const DiffusionTensorImagePointer & tensorImage,
      int term,
      const SpacingType & spacing,
      const typename OutputImageType::SizeType & radius );

  /** Allocate and populate the images of multiplication vectors that the
   *  div( T \grad( u ) ) values are multiplied by.  Allocate and populate
   *  all or
   *  some of the multiplication vector images in derived classes.
   *  Otherwise,
   *  default to e_l, where e_l is the lth canonical unit vector. */
  virtual void ComputeMultiplicationVectorImages( void ) {}

  /** Initialize the state of the filter before each iteration. */
  virtual void InitializeIteration( void ) override;

  /** Updates the deformation vector component images on each iteration. */
  virtual void UpdateDeformationComponentImages( OutputImageType * output );

  /** Computes the first- and second-order partial derivatives of the
   *  deformation component images on each iteration.  Override in derived
   *  classes if the deformation components image pointers are not unique,
   *  to
   *  avoid computing the same derivatives multiple times. */
  virtual void ComputeDeformationComponentDerivativeImages( void );

  /** Helper to compute the first- and second-order partial derivatives of
   * the
   *  deformation component images, using the
   *  ThreadedComputeDeformationComponentDerivativeImageHelper() method and
   *  a
   *  multithreading mechanism.
   *  \sa ThreadedComputeDeformationComponentDerivativeImageHelper */
  virtual void ComputeDeformationComponentDerivativeImageHelper(
      DeformationVectorComponentImagePointer & deformationComponentImage,
      int term,
      int dimension,
      SpacingType & spacing,
      typename OutputImageType::SizeType & radius );

  /** Does the actual work of computing the first- and second-order partial
   *  derivatives of the deformation component images over regions supplied
   *  by the multithreading mechanism.
   *  \sa ComputeDeformationComponentDerivativeImageHelper
   *  \sa ComputeDeformationComponentDerivativeImageHelperThreadedCallback
   *  */
   virtual void ThreadedComputeDeformationComponentDerivativeImageHelper(
       const DeformationVectorComponentImagePointer
         & deformationComponentImage,
       const ThreadDeformationVectorComponentImageRegionType
         & deformationVectorComponenntRegionToProcess,
       const ThreadScalarDerivativeImageRegionType
         & scalarDerivativeRegionToProcess,
       const ThreadTensorDerivativeImageRegionType
         & tensorDerivativeRegionToProcess,
       int term,
       int dimension,
       const SpacingType & spacing,
       const typename OutputImageType::SizeType & radius ) const;

  /** Get a diffusion tensor image */
  DiffusionTensorImageType * GetDiffusionTensorImage( int index ) const
    {
    assert( index < this->GetNumberOfTerms() );
    return this->m_DiffusionTensorImages[index];
    }

  /** Get an image of the diffusion tensor derivatives */
  TensorDerivativeImageType * GetDiffusionTensorDerivativeImage( int index )
      const
    {
    assert( index < this->GetNumberOfTerms() );
    return this->m_DiffusionTensorDerivativeImages[index];
    }

  /** Set/get an image of the deformation field components */
  DeformationFieldType * GetDeformationComponentImage( int index ) const
    {
    assert( index < this->GetNumberOfTerms() );
    return this->m_DeformationComponentImages[index];
    }
  void SetDeformationComponentImage( int index, DeformationFieldType *
    comp )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( comp );
    this->m_DeformationComponentImages[index] = comp;
    }

  /** Set/get a multiplication vectors image. */
  void SetMultiplicationVectorImage( int index,
    int dimension,
    DeformationFieldType * mult )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( dimension < static_cast< int >( ImageDimension ) );
    this->m_MultiplicationVectorImageArrays[index][dimension] = mult;
    }
  DeformationFieldType * GetMultiplicationVectorImage( int index,
    int dimension )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( dimension < static_cast< int >( ImageDimension ) );
    return this->m_MultiplicationVectorImageArrays[index][dimension];
    }

  /** Set/Get a first-order deformation component derivative. */
  void SetDeformationComponentFirstOrderDerivative(
      int index,
      int dimension,
      ScalarDerivativeImageType * deriv )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( dimension < static_cast< int >( ImageDimension ) );
    this->m_DeformationComponentFirstOrderDerivativeArrays[index][dimension]
        = deriv;
    }
  ScalarDerivativeImageType * GetDeformationComponentFirstOrderDerivative(
      int index,
      int dimension )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( dimension < static_cast< int >( ImageDimension ) );
    return this->m_DeformationComponentFirstOrderDerivativeArrays
        [index][dimension];
    }

  /** Set/Get a second-order deformation component derivative. */
  void SetDeformationComponentSecondOrderDerivative(
      int index,
      int dimension,
      TensorDerivativeImageType * deriv )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( dimension < static_cast< int >( ImageDimension ) );
    this->m_DeformationComponentSecondOrderDerivativeArrays[index][dimension]
        = deriv;
    }
  TensorDerivativeImageType * GetDeformationComponentSecondOrderDerivative(
      int index,
      int dimension )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( dimension < static_cast< int >( ImageDimension ) );
    return this->m_DeformationComponentSecondOrderDerivativeArrays
        [index][dimension];
    }

  struct UpdateMetricsIntermediateStruct
    {
    int     NumberOfPixelsProcessed;
    double  SumOfSquaredTotalUpdateMagnitude;
    double  SumOfSquaredIntensityDistanceUpdateMagnitude;
    double  SumOfSquaredRegularizationUpdateMagnitude;
    double  SumOfTotalUpdateMagnitude;
    double  SumOfIntensityDistanceUpdateMagnitude;
    double  SumOfRegularizationUpdateMagnitude;

    void zero( void )
      {
      NumberOfPixelsProcessed = 0;
      SumOfSquaredTotalUpdateMagnitude = 0.0;
      SumOfSquaredIntensityDistanceUpdateMagnitude = 0.0;
      SumOfSquaredRegularizationUpdateMagnitude = 0.0;
      SumOfTotalUpdateMagnitude = 0.0;
      SumOfIntensityDistanceUpdateMagnitude = 0.0;
      SumOfRegularizationUpdateMagnitude = 0.0;
      }

    void copyFrom( const UpdateMetricsIntermediateStruct & rhs )
      {
      NumberOfPixelsProcessed = rhs.NumberOfPixelsProcessed;
      SumOfSquaredTotalUpdateMagnitude =
        rhs.SumOfSquaredTotalUpdateMagnitude;
      SumOfSquaredIntensityDistanceUpdateMagnitude
          = rhs.SumOfSquaredIntensityDistanceUpdateMagnitude;
      SumOfSquaredRegularizationUpdateMagnitude
          = rhs.SumOfSquaredRegularizationUpdateMagnitude;
      SumOfTotalUpdateMagnitude = rhs.SumOfTotalUpdateMagnitude;
      SumOfIntensityDistanceUpdateMagnitude
          = rhs.SumOfIntensityDistanceUpdateMagnitude;
      SumOfRegularizationUpdateMagnitude
          = rhs.SumOfRegularizationUpdateMagnitude;
      }

    void difference( const UpdateMetricsIntermediateStruct & lhs,
                     const UpdateMetricsIntermediateStruct & rhs )
      {
      SumOfSquaredTotalUpdateMagnitude =
        lhs.SumOfSquaredTotalUpdateMagnitude
          - rhs.SumOfSquaredTotalUpdateMagnitude;
      SumOfSquaredIntensityDistanceUpdateMagnitude =
        lhs.SumOfSquaredIntensityDistanceUpdateMagnitude
          - rhs.SumOfSquaredIntensityDistanceUpdateMagnitude;
      SumOfSquaredRegularizationUpdateMagnitude =
        lhs.SumOfSquaredRegularizationUpdateMagnitude
          - rhs.SumOfSquaredRegularizationUpdateMagnitude;
      SumOfTotalUpdateMagnitude =
        lhs.SumOfTotalUpdateMagnitude - rhs.SumOfTotalUpdateMagnitude;
      SumOfIntensityDistanceUpdateMagnitude =
        lhs.SumOfIntensityDistanceUpdateMagnitude
          - rhs.SumOfIntensityDistanceUpdateMagnitude;
      SumOfRegularizationUpdateMagnitude =
        lhs.SumOfRegularizationUpdateMagnitude
          - rhs.SumOfRegularizationUpdateMagnitude;
      }
    }; // End struct UpdateMetricsIntermediateStruct

  struct UpdateMetricsStruct
    {
    UpdateMetricsIntermediateStruct IntermediateStruct;
    double                          RMSTotalUpdateMagnitude;
    double                          RMSIntensityDistanceUpdateMagnitude;
    double                          RMSRegularizationUpdateMagnitude;
    double                          MeanTotalUpdateMagnitude;
    double                          MeanIntensityDistanceUpdateMagnitude;
    double                          MeanRegularizationUpdateMagnitude;

    void zero( void )
      {
      IntermediateStruct.zero();
      RMSTotalUpdateMagnitude = 0.0;
      RMSIntensityDistanceUpdateMagnitude = 0.0;
      RMSRegularizationUpdateMagnitude = 0.0;
      MeanTotalUpdateMagnitude = 0.0;
      MeanIntensityDistanceUpdateMagnitude = 0.0;
      MeanRegularizationUpdateMagnitude = 0.0;
      }

    void copyFrom( const UpdateMetricsStruct & rhs )
      {
      IntermediateStruct.copyFrom( rhs.IntermediateStruct );
      RMSTotalUpdateMagnitude = rhs.RMSTotalUpdateMagnitude;
      RMSIntensityDistanceUpdateMagnitude
          = rhs.RMSIntensityDistanceUpdateMagnitude;
      RMSRegularizationUpdateMagnitude =
        rhs.RMSRegularizationUpdateMagnitude;
      MeanTotalUpdateMagnitude = rhs.MeanTotalUpdateMagnitude;
      MeanIntensityDistanceUpdateMagnitude
          = rhs.MeanIntensityDistanceUpdateMagnitude;
      MeanRegularizationUpdateMagnitude =
        rhs.MeanRegularizationUpdateMagnitude;
      }

    // Does not calculate for intermediate struct
    void difference( const UpdateMetricsStruct & lhs,
                     const UpdateMetricsStruct & rhs )
      {
      RMSTotalUpdateMagnitude =
        lhs.RMSTotalUpdateMagnitude - rhs.RMSTotalUpdateMagnitude;
      RMSIntensityDistanceUpdateMagnitude =
        lhs.RMSIntensityDistanceUpdateMagnitude
          - rhs.RMSIntensityDistanceUpdateMagnitude;
      RMSRegularizationUpdateMagnitude =
        lhs.RMSRegularizationUpdateMagnitude
          - rhs.RMSRegularizationUpdateMagnitude;
      MeanTotalUpdateMagnitude =
        lhs.MeanTotalUpdateMagnitude - rhs.MeanTotalUpdateMagnitude;
      MeanIntensityDistanceUpdateMagnitude =
        lhs.MeanIntensityDistanceUpdateMagnitude
          - rhs.MeanIntensityDistanceUpdateMagnitude;
      MeanRegularizationUpdateMagnitude =
        lhs.MeanRegularizationUpdateMagnitude
          - rhs.MeanRegularizationUpdateMagnitude;
      }
    }; // End strut UpdateMetricsStruct

  /** This method populates an update buffer with changes for each pixel
   * in the
   * output.  Uses CalculateChangeGradient to determine
   * the gradient displacement field.
   * The update buffer does not include the global scaling
   * parameter, which is incorporated in ApplyUpdate
   * Return value is a time step to be used for the update.
   * \sa CalculateChangeGradient */
  virtual TimeStepType CalculateChange( void ) override;

  /** Inherited from superclass - do not call this function! */
  TimeStepType ThreadedCalculateChange(
      const ThreadRegionType & regionToProcess, ThreadIdType threadId ) override;

  /** This method populates an update buffer with changes for each pixel
   * in the
   * output when computing the gradient, using the
   * ThreadedCalculateChange() method and a multithreading
   * mechanism. Return value is a time step to be used for the update.
   * \sa ThreadedCalculateChangeGradient */
  virtual TimeStepType CalculateChangeGradient( void );

  /** Does the actual work of calculating the gradient part of the line
   * search
   * over a region supplied by the multithreading mechanism.
   * \sa CalculateChangeGradient
   * \sa CalculateChangeGradientThreaderCallback */
  virtual TimeStepType ThreadedCalculateChangeGradient(
      const ThreadRegionType & regionToProcess,
      const ThreadDiffusionTensorImageRegionType & tensorRegionToProcess,
      const ThreadTensorDerivativeImageRegionType
        & tensorDerivativeRegionToProcess,
      const ThreadScalarDerivativeImageRegionType
        & scalarDerivativeRegionToProcess,
      const ThreadStoppingCriterionMaskImageRegionType
        & stoppingCriterionMaskRegionToProcess,
      UpdateMetricsIntermediateStruct & updateMetricsIntermediate,
      int threadId );

  /** Calculates the total, intensity distance and regularization energies,
   *  using the ThreadedCalculateEnergies() method and a multithreading
   *  mechanism.  The stepSize parameter is a uniform scaling parameter
   *  with which to scale the displacement field for use with the line
   *  search.
   * \sa ThreadedCalculateEnergies */
  virtual void CalculateEnergies( EnergiesStruct & energies,
                                  OutputImageType * outputField );

  /** Updates the intermediate update statistics ( sum-of-squared and sum-of
   *  statistics, incorporating the time step, and
   *  computes the RMS and mean update statistics */
  virtual void UpdateUpdateStatistics( TimeStepType stepSize );

  /** Does the actual work of calculating the intensity distance and
   *  regularization energies over a region supplied by the multithreading
   * mechanism.
   * \sa CalculateEnergies
   * \sa CalculateEnergiesThreaderCallback */
  virtual void ThreadedCalculateEnergies(
    const OutputImagePointer & output,
    const ThreadRegionType & regionToProcess,
    const ThreadDiffusionTensorImageRegionType & tensorRegionToProcess,
    const ThreadScalarDerivativeImageRegionType &
      scalarDerivativeRegionToProcess,
    const ThreadStoppingCriterionMaskImageRegionType &
      stoppingCriterionMaskRegionToProcess,
    double & intensityDistanceEnergy,
    double & regularizationEnergy,
    int threadId );

  /** This method applies changes from the update buffer to the output,
   * using
   * the ThreadedApplyUpdate() method and a multithreading mechanism.
   * "dt" is
   * the time step to use for the update of each pixel.  Also multiplies
   * each
   * voxel in the update buffer by the scaling value from the line search.
   * \sa ThreadedApplyUpdate */
  virtual void ApplyUpdate( const TimeStepType & dt ) override;
  virtual void ApplyUpdate( TimeStepType dt, OutputImagePointer
    outputImage );

  /** Inherited from superclass - do not call this function! */
  virtual void ThreadedApplyUpdate( const TimeStepType & dt,
    const ThreadRegionType & regionToProcess,
    ThreadIdType threadId ) override;

  /**  Does the actual work of updating the output from the UpdateContainer
   * over an output region supplied by the multithreading mechanism.
   *  \sa ApplyUpdate
   *  \sa ApplyUpdateThreaderCallback */
  virtual void ThreadedApplyUpdate( OutputImagePointer & outputImage,
      TimeStepType dt, const ThreadRegionType & regionToProcess,
      ThreadIdType threadId );

  /** Create the registration function, with default parameters for
    * ComputeRegularizationTerm and ComputeIntensityDistanceTerm. */
  virtual void CreateRegistrationFunction( void );

  /** Get the registration function pointer */
  virtual RegistrationFunctionType * GetRegistrationFunctionPointer(
    void ) const;

  /** Allocate the update buffer. */
  virtual void AllocateUpdateBuffer( void ) override;

  /** Get the update buffer. */
  virtual UpdateBufferType * GetUpdateBuffer( void ) override
    { return m_UpdateBuffer; }

  /** This method is called after ApplyUpdate() to print out energy and RMS
   * change metrics and evaluate the stopping conditions. */
  virtual void PostProcessIteration( TimeStepType stepSize );

private:
  // Purposely not implemented
  DiffusiveRegistrationFilter( const Self& );
  void operator=( const Self& ); // Purposely not implemented

  /** Structure for passing information into static callback methods.  Used
   * in
   * the subclasses threading mechanisms. */
  struct DenseFDThreadStruct
    {
    DiffusiveRegistrationFilter *Filter;
    OutputImagePointer OutputImage;
    TimeStepType TimeStep;
    TimeStepType *TimeStepList;
    itk::BooleanStdVectorType *ValidTimeStepList;
    }; // End struct DenseFDThreadStruct

  /** Structure for passing information into static callback methods.  Used
   * in
   *  the threading mechanism for
   *  ComputeDeformationComponentDerivativeImageHelper. */
  struct ComputeDeformationComponentDerivativeImageHelperThreadStruct
    {
    DiffusiveRegistrationFilter * Filter;
    DeformationVectorComponentImagePointer DeformationComponentImage;
    int Term;
    int Dimension;
    SpacingType Spacing;
    typename OutputImageType::SizeType Radius;
    };
  // End struct ComputeDeformationComponentDerivativeImageHelperThreadStruct

  /** Structure for passing information into static callback methods.  Used
   * in
   *  the threading mechanism for CalculateChangeGradient. */
  struct CalculateChangeGradientThreadStruct
    {
    DiffusiveRegistrationFilter *Filter;
    TimeStepType TimeStep;
    std::vector< TimeStepType > TimeStepList;
    itk::BooleanStdVectorType ValidTimeStepList;
    UpdateMetricsIntermediateStruct *UpdateMetricsIntermediate;
    }; // End struct CalculateChangeGradientThreadStruct

  /** Structure for passing information into static callback methods.
   * Used in
   *  the threading mechanism for CalculateEnergies. */
  struct CalculateEnergiesThreadStruct
    {
    DiffusiveRegistrationFilter * Filter;
    OutputImagePointer OutputImage;
    double * IntensityDistanceEnergies;
    double * RegularizationEnergies;
    }; // End struct CalculateEnergiesThreadStruct

  /** This callback method uses ImageSource::SplitRequestedRegion to
   * acquire an
   * output region that it passes to ThreadedApplyUpdate for processing. */
  static ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
    ApplyUpdateThreaderCallback( void *arg );

  /** This callback method uses SplitUpdateContainer to acquire a region
   * which it then passes to ThreadedCalculateChange for processing. */
  static ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
    CalculateChangeGradientThreaderCallback( void *arg );

  /** This callback method uses SplitUpdateContainer to acquire a region
   * which it then passes to
   * ThreadedComputeDeformationComponentDerivativeImageHelper
   * for processing. */
  static ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
   ComputeDeformationComponentDerivativeImageHelperThreaderCallback(
     void *arg );

  /** This callback method uses SplitUpdateContainer to acquire a region
   * which it then passes to ThreadedComputeEnergies for processing. */
  static ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
    CalculateEnergiesThreaderCallback( void *arg );

  TimeStepType                              m_OriginalTimeStep;

  /** The buffer that holds the updates for an iteration of algorithm,
   * without
   *  the scaling provided by the line search. */
  typename UpdateBufferType::Pointer        m_UpdateBuffer;

  /** Images storing information we will need for each voxel on every
   *  registration iteration */
  DiffusionTensorImageArrayType       m_DiffusionTensorImages;
  TensorDerivativeImageVectorType     m_DiffusionTensorDerivativeImages;
  DeformationFieldArrayType           m_DeformationComponentImages;

  ScalarDerivativeImageArrayVectorType
    m_DeformationComponentFirstOrderDerivativeArrays;
  TensorDerivativeImageArrayVectorType
    m_DeformationComponentSecondOrderDerivativeArrays;
  DeformationVectorImageArrayVectorType
    m_MultiplicationVectorImageArrays;

  /** Variables for multiresolution registration.  Current level can be
   * detected
   *  as Initialize() is called on each new level. */
  unsigned int                              m_CurrentLevel;

  /** Relative weightings between the intensity distance and regularization
   *  update terms.  Stored in a vector so that the user can provide
   *  different
   *  weightings per multiresolution level. */
  std::vector< double >                     m_RegularizationWeightings;

  /** Template used to calculate member images */
  FixedImagePointer                         m_HighResolutionTemplate;

  /** Mask on the stopping criterion computation */
  StoppingCriterionMaskPointer   m_HighResolutionStoppingCriterionMask;
  StoppingCriterionMaskPointer   m_StoppingCriterionMask;

  /** Parameters for stopping criterion */
  unsigned int                   m_StoppingCriterionEvaluationPeriod;
  double                         m_StoppingCriterionMaxTotalEnergyChange;

  /** Parameters for energies and update magnitude metrics */
  EnergiesStruct                            m_Energies;
  EnergiesStruct                            m_PreviousEnergies;
  UpdateMetricsStruct                       m_UpdateMetrics;
  UpdateMetricsStruct                       m_PreviousUpdateMetrics;

}; // End class DiffusiveRegistrationFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeDiffusiveRegistrationFilter.hxx"
#endif

#endif // End !defined( __itktubeDiffusiveRegistrationFilter_h )
