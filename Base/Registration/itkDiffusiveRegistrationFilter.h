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
#ifndef __itkDiffusiveRegistrationFilter_h
#define __itkDiffusiveRegistrationFilter_h

#include "itkPDEDeformableRegistrationFilter.h"
#include "itkAnisotropicDiffusiveRegistrationFunction.h"

#include "itkVectorIndexSelectionCastImageFilter.h"

namespace itk
{

/** \class itkDiffusiveRegistrationFilter
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
 * example, these regularizers can accomodate deformation field discontinuities
 * that are expected when considering sliding motion.
 *
 * The update term is composed of two parts: an intensity distance term and a
 * regularization term.  The intensity distance term uses the sum of square
 * difference metric, so this registration algorithm is appropriate for
 * monomodal image registration term only.
 *
 * The update term for the regularization will be of the form:
 * div(T1*\grad(u1))v1 + div(T2*\grad(u2))v2 + ... + div(TN*\grad(uN))vN
 * where the types are:
 * - T1..TN are diffusion tensors
 * - u1..uN are deformation vectors
 * - v1..vN are deformation vectors
 * It is assumed that T's and v's are constant throughout the registration, and
 * so can be precomputed once, while u's must be updated on each registration
 * iteration.
 *
 * This base class implements the diffusive regularization.  Algorithms
 * implementing anisotropic regularization should derive it and override the
 * following functions:
 * - GetNumberOfTerms(): returns the number of div(T*\grad(u))v terms
 * - ComputeDiffusionTensorImages(): allocate and populate the T images
 * - InitializeDeformationComponentAndDerivativeImages(): allocate the u images
 * and their derivatives
 * - ComputeMultiplicationVectorImages(): allocate and populate the v images
 * - UpdateDeformationComponentImages(): update the u images at each iteration
 * See itkAnisotropicDiffusiveRegistrationFilter for an example derived filter.
 *
 * See: D.F. Pace et al., Deformable image registration of sliding organs using
 * anisotropic diffusive regularization, ISBI 2011.
 *
 * This class is templated over the type of the fixed image, the type of the
 * moving image and the type of the deformation field.
 *
 * \sa itkAnisotropicDiffusiveRegistrationFunction
 * \sa itkAnisotropicDiffusiveRegistrationFilter
 * \ingroup DeformableImageRegistration
 * \ingroup MultiThreaded
 */

template < class TFixedImage, class TMovingImage, class TDeformationField >
class ITK_EXPORT DiffusiveRegistrationFilter
  : public PDEDeformableRegistrationFilter< TFixedImage,
                                           TMovingImage,
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
    * itkNewMacro(), but the type of the registration function depends on the
    * type of this object.  Can't call the overridden function
    * CreateRegistrationFunction() from the base class constructor, so we'll
    * call it here. Derived classes should use this instead of itkNewMacro().*/
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, PDEDeformableRegistrationFilter);

  /** Inherit some parameters from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Convenient typedefs from the superclass. */
  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::FixedImagePointer        FixedImagePointer;
  typedef typename Superclass::MovingImageType          MovingImageType;
  typedef typename Superclass::MovingImagePointer       MovingImagePointer;
  typedef typename Superclass::DeformationFieldType     DeformationFieldType;
  typedef typename Superclass::DeformationFieldPointer  DeformationFieldPointer;
  typedef typename Superclass::TimeStepType             TimeStepType;
  typedef typename Superclass::FiniteDifferenceFunctionType
      FiniteDifferenceFunctionType;

  /** Typedefs used in multithreading */
  typedef typename Superclass::OutputImageType          OutputImageType;
  typedef typename Superclass::OutputImagePointer       OutputImagePointer;
  typedef typename Superclass::UpdateBufferType         UpdateBufferType;
  typedef typename UpdateBufferType::RegionType         ThreadRegionType;

  /** Output image and update buffer types */
  typedef itk::ImageRegionIterator< OutputImageType >   OutputImageRegionType;
  typedef typename FiniteDifferenceFunctionType::NeighborhoodType
      NeighborhoodType;
  typedef itk::ImageRegionIterator< UpdateBufferType >  UpdateBufferRegionType;

  /** The registration function type */
  typedef AnisotropicDiffusiveRegistrationFunction
      < FixedImageType, MovingImageType, DeformationFieldType >
      RegistrationFunctionType;
  typedef typename RegistrationFunctionType::RegularizationFunctionType
      RegularizationFunctionType;
  typedef typename RegistrationFunctionType::RegularizationFunctionPointer
      RegularizationFunctionPointer;
  typedef typename RegistrationFunctionType::SpacingType SpacingType;

  /** Deformation component types (i.e. component of a deformation field,
   *  still a vector */
  typedef std::vector< DeformationFieldPointer >
      DeformationFieldArrayType;
  typedef typename RegistrationFunctionType::DeformationVectorType
      DeformationVectorType;

  /** Deformation vector component types (i.e. scalar within a vector) */
  typedef typename RegistrationFunctionType::DeformationVectorComponentType
      DeformationVectorComponentType;
  typedef typename RegistrationFunctionType::DeformationVectorComponentImageType
      DeformationVectorComponentImageType;
  typedef typename DeformationVectorComponentImageType::Pointer
      DeformationVectorComponentImagePointer;
  typedef typename
      itk::FixedArray< DeformationVectorComponentImagePointer, ImageDimension >
      DeformationComponentImageArrayType;
  typedef typename
      RegistrationFunctionType::DeformationVectorComponentNeighborhoodType
      DeformationVectorComponentNeighborhoodType;

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
  typedef typename RegistrationFunctionType::DeformationVectorImageRegionType
      DeformationVectorImageRegionType;
  typedef typename
      RegistrationFunctionType::DeformationVectorImageRegionArrayVectorType
      DeformationVectorImageRegionArrayVectorType;

  /** Types for vector component extractor */
  typedef itk::VectorIndexSelectionCastImageFilter
      < DeformationFieldType, DeformationVectorComponentImageType >
      VectorIndexSelectionFilterType;

  /** Convenience functions to set/get the registration functions timestep. */
  void SetTimeStep( const TimeStepType & t )
    { this->GetRegistrationFunctionPointer()->SetTimeStep( t ); }
  const TimeStepType& GetTimeStep() const
    { return this->GetRegistrationFunctionPointer()->GetTimeStep(); }

  /** Set/get whether to compute the motion field regularization term
   *  Default: true */
  void SetComputeRegularizationTerm( bool compute )
    { this->GetRegistrationFunctionPointer()->
      SetComputeRegularizationTerm( compute ); }
  bool GetComputeRegularizationTerm() const
    { return this->GetRegistrationFunctionPointer()->
      GetComputeRegularizationTerm(); }

  /** Set/get whether to compute the intensity distance term
   *  Default: true */
  void SetComputeIntensityDistanceTerm( bool compute )
    { this->GetRegistrationFunctionPointer()->
      SetComputeIntensityDistanceTerm( compute ); }
  bool GetComputeIntensityDistanceTerm() const
    { return this->GetRegistrationFunctionPointer()->
      GetComputeIntensityDistanceTerm(); }

  /** The number of div(T\grad(u))v terms we sum for the regularizer.
   *  Reimplement in derived classes. */
  virtual int GetNumberOfTerms() const
    { return 1; }

  /** Set/get a pointer to an image that is to be used for the template when
   *  computing member images.  This is usually the original fixed image.  The
   *  attributes of this filter's output are used if the high resolution
   *  template is not set.  For proper behavior, you must set this if using a
   *  multiresolution registration. */
  virtual void SetHighResolutionTemplate( FixedImageType * templateImage )
    { m_HighResolutionTemplate = templateImage; }
  virtual FixedImageType * GetHighResolutionTemplate()
    { return m_HighResolutionTemplate; }

protected:
  DiffusiveRegistrationFilter();
  virtual ~DiffusiveRegistrationFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Handy for array indexing. */
  enum DivTerm { GAUSSIAN };

  /** Initialization occuring before the registration iterations begin. */
  virtual void Initialize();

  /** Allocate images used during the registration. */
  virtual void AllocateImageMembers();

  /** Allocate the deformation component images and their derivative images.
   *  (which may be updated throughout the registration). Reimplement in derived
   *  classes. */
  virtual void InitializeDeformationComponentAndDerivativeImages();

  /** Allocate and populate the diffusion tensor images.
   *  Reimplement in derived classes. */
  virtual void ComputeDiffusionTensorImages();

  /** Computes the first-order partial derivatives of the diffusion tensor
   *  images.  Override in derived classes if the diffusion tensor image
   *  pointers are not unique, to avoid computing the derivatives multiple
   *  times. */
  virtual void ComputeDiffusionTensorDerivativeImages();

  /** Helper to compute the first-order partial derivatives of the diffusion
   *  tensor images */
  virtual void ComputeDiffusionTensorDerivativeImageHelper(
      const DiffusionTensorImagePointer & tensorImage,
      int term,
      const SpacingType & spacing,
      const typename OutputImageType::SizeType & radius ) const;

  /** Allocate and populate the images of multiplication vectors that the
   *  div(T \grad(u)) values are multiplied by.  Allocate and populate all or
   *  some of the multiplication vector images in derived classes.  Otherwise,
   *  default to e_l, where e_l is the lth canonical unit vector. */
  virtual void ComputeMultiplicationVectorImages() {};

  /** Initialize the state of the filter before each iteration. */
  virtual void InitializeIteration();

  /** Updates the deformation vector component images on each iteration. */
  virtual void UpdateDeformationComponentImages() {};

  /** Computes the first- and second-order partial derivatives of the
   *  deformation component images on each iteration.  Override in derived
   *  classes if the deformation components image pointers are not unique, to
   *  avoid computing the same derivatives multiple times. */
  virtual void ComputeDeformationComponentDerivativeImages();

  /** Helper to compute the first- and second-order partial derivatives of the
   *  deformation component images */
  virtual void ComputeDeformationComponentDerivativeImageHelper(
      const DeformationVectorComponentImagePointer & deformationComponentImage,
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
  void SetDeformationComponentImage( int index, DeformationFieldType * comp )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( comp );
    this->m_DeformationComponentImages[index] = comp;
    }

  /** Set/get a multiplication vectors image. */
  void SetMultiplicationVectorImage(
      int index, int dimension, DeformationFieldType * mult )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( dimension < ImageDimension );
    this->m_MultiplicationVectorImageArrays[index][dimension] = mult;
    }
  DeformationFieldType * GetMultiplicationVectorImage( int index,
                                                             int dimension )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( dimension < ImageDimension );
    return this->m_MultiplicationVectorImageArrays[index][dimension];
    }

  /** Set/Get a first-order deformation component derivative. */
  void SetDeformationComponentFirstOrderDerivative(
      int index, int dimension, ScalarDerivativeImageType * deriv )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( dimension < ImageDimension );
    this->m_DeformationComponentFirstOrderDerivativeArrays[index][dimension]
        = deriv;
    }
  ScalarDerivativeImageType * GetDeformationComponentFirstOrderDerivative(
      int index, int dimension )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( dimension < ImageDimension );
    return this->m_DeformationComponentFirstOrderDerivativeArrays
        [index][dimension];
    }

  /** Set/Get a second-order deformation component derivative. */
  void SetDeformationComponentSecondOrderDerivative(
      int index, int dimension, TensorDerivativeImageType * deriv )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( dimension < ImageDimension );
    this->m_DeformationComponentSecondOrderDerivativeArrays[index][dimension]
        = deriv;
    }
  TensorDerivativeImageType * GetDeformationComponentSecondOrderDerivative(
      int index, int dimension )
    {
    assert( index < this->GetNumberOfTerms() );
    assert( dimension < ImageDimension );
    return this->m_DeformationComponentSecondOrderDerivativeArrays
        [index][dimension];
    }

  /** Extracts the x, y, z components of a deformation field. */
  void ExtractXYZComponentsFromDeformationField(
      const OutputImageType * deformationField,
      DeformationComponentImageArrayType & deformationComponentImages ) const;

  /** This method populates an update buffer with changes for each pixel in the
   * output, using the ThreadedCalculateChange() method and a multithreading
   * mechanism. Return value is a time step to be used for the update.
   * \sa ThreadedCalculateChange */
  virtual TimeStepType CalculateChange();

  /** Inherited from superclass - do not call this function!  Call the other
   *  ThreadedCalculateChange function instead */
  TimeStepType ThreadedCalculateChange(
      const ThreadRegionType & regionToProcess, int threadId );

  /** Does the actual work of calculating change over a region supplied by
   * the multithreading mechanism.
   * \sa CalculateChange
   * \sa CalculateChangeThreaderCallback */
  virtual TimeStepType ThreadedCalculateChange(
      const ThreadRegionType & regionToProcess,
      const ThreadDiffusionTensorImageRegionType & tensorRegionToProcess,
      const ThreadTensorDerivativeImageRegionType
        & tensorDerivativeRegionToProcess,
      const ThreadScalarDerivativeImageRegionType
        & scalarDerivativeRegionToProcess,
      int threadId );

  /** This method applies changes from the update buffer to the output, using
   * the ThreadedApplyUpdate() method and a multithreading mechanism.  "dt" is
   * the time step to use for the update of each pixel.
   * \sa ThreadedApplyUpdate */
  virtual void ApplyUpdate( TimeStepType dt );

  /**  Does the actual work of updating the output from the UpdateContainer over
   *  an output region supplied by the multithreading mechanism.
   *  \sa ApplyUpdate
   *  \sa ApplyUpdateThreaderCallback */
  virtual void ThreadedApplyUpdate( TimeStepType dt,
                                    const ThreadRegionType & regionToProcess,
                                    int threadId );

  /** Create the registration function, with default parameters for
    * ComputeRegularizationTerm and ComputeIntensityDistanceTerm. */
  virtual void CreateRegistrationFunction();

  /** Get the registration function pointer */
  virtual RegistrationFunctionType * GetRegistrationFunctionPointer() const;

  /** Allocate the update buffer. */
  virtual void AllocateUpdateBuffer();

  /** Get the update buffer. */
  virtual UpdateBufferType * GetUpdateBuffer()
    { return m_UpdateBuffer; }

  /** Helper function to allocate an image based on a template */
  template< class UnallocatedImagePointer, class TemplateImagePointer >
  void AllocateSpaceForImage( UnallocatedImagePointer & image,
                              const TemplateImagePointer & templateImage )
  const;

  /** Helper function to check whether the attributes of an image match a
    * template */
  template< class CheckedImageType, class TemplateImageType >
  bool CompareImageAttributes( const CheckedImageType * image,
                               const TemplateImageType * templateImage )
  const;

  /** Resamples an image to a template using nearest neighbor interpolation */
  template< class ResampleImageType, class TemplateImageType >
  void ResampleImageNearestNeighbor(
      const ResampleImageType * highResolutionImage,
      const TemplateImageType * templateImage,
      ResampleImageType * resampledImage ) const;

  /** Resamples an image to a template using linear interpolation */
  template< class ResampleImageType, class TemplateImageType >
  void ResampleImageLinear(
      const ResampleImageType * highResolutionImage,
      const TemplateImageType * templateImage,
      ResampleImageType * resampledImage ) const;

  /** Resamples a vector image to a template using linear interpolation */
  template< class VectorResampleImageType, class TemplateImageType >
  void VectorResampleImageLinear(
      const VectorResampleImageType * highResolutionImage,
      const TemplateImageType * templateImage,
      VectorResampleImageType * resampledImage ) const;

private:
  // Purposely not implemented
  DiffusiveRegistrationFilter(const Self&);
  void operator=(const Self&); // Purposely not implemented

  /** Structure for passing information into static callback methods.  Used in
   * the subclasses threading mechanisms. */
  struct DenseFDThreadStruct
    {
    DiffusiveRegistrationFilter *Filter;
    TimeStepType TimeStep;
    TimeStepType *TimeStepList;
    bool *ValidTimeStepList;
    };

  /** This callback method uses ImageSource::SplitRequestedRegion to acquire an
   * output region that it passes to ThreadedApplyUpdate for processing. */
  static ITK_THREAD_RETURN_TYPE ApplyUpdateThreaderCallback( void *arg );

  /** This callback method uses SplitUpdateContainer to acquire a region
  * which it then passes to ThreadedCalculateChange for processing. */
  static ITK_THREAD_RETURN_TYPE CalculateChangeThreaderCallback( void *arg );

  /** The buffer that holds the updates for an iteration of algorithm. */
  typename UpdateBufferType::Pointer  m_UpdateBuffer;

  /** Images storing information we will need for each voxel on every
   *  registration iteration */
  DiffusionTensorImageArrayType             m_DiffusionTensorImages;
  TensorDerivativeImageVectorType           m_DiffusionTensorDerivativeImages;
  DeformationFieldArrayType                 m_DeformationComponentImages;
  ScalarDerivativeImageArrayVectorType
      m_DeformationComponentFirstOrderDerivativeArrays;
  TensorDerivativeImageArrayVectorType
      m_DeformationComponentSecondOrderDerivativeArrays;
  DeformationVectorImageArrayVectorType     m_MultiplicationVectorImageArrays;

  /** Template used to calculate member images */
  FixedImagePointer                         m_HighResolutionTemplate;
};

/** Struct to simply get the face list and an iterator over the face list
 *  when processing an image.  Designed for use with SmartPointers. */
template< class ImageType >
struct FaceStruct
  {
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
      < typename ImageType::ObjectType > FaceCalculatorType;
  typedef typename FaceCalculatorType::FaceListType FaceListType;
  typedef typename FaceListType::iterator FaceListIteratorType;

  FaceStruct()
    {
    numberOfTerms = 0;
    }

  FaceStruct( const ImageType& image,
              typename ImageType::ObjectType::RegionType region,
              typename ImageType::ObjectType::SizeType radius )
    {
    numberOfTerms = 0;
    if( image.GetPointer() )
      {
      faceLists.push_back( faceCalculator( image, region, radius ) );
      numberOfTerms = 1;
      }
    }

  FaceStruct( const std::vector< ImageType >& images,
              typename ImageType::ObjectType::RegionType region,
              typename ImageType::ObjectType::SizeType radius )
    {
    numberOfTerms = 0;
    for( int i = 0; i < (int) images.size(); i++ )
      {
      if( images[i].GetPointer() )
        {
        faceLists.push_back( faceCalculator( images[i], region, radius ) );
        numberOfTerms++;
        }
      }
    }

  template< unsigned int VLength >
  FaceStruct( const itk::FixedArray< ImageType, VLength >& images,
              typename ImageType::ObjectType::RegionType region,
              typename ImageType::ObjectType::SizeType radius )
    {
    numberOfTerms = 0;
    for( int i = 0; i < (int) images.Size(); i++ )
      {
      if( images[i].GetPointer() )
        {
        faceLists.push_back( faceCalculator( images[i], region, radius ) );
        numberOfTerms++;
        }
      }
    }

  template< unsigned int VLength >
  FaceStruct( const
              std::vector< itk::FixedArray< ImageType, VLength > > &images,
              typename ImageType::ObjectType::RegionType region,
              typename ImageType::ObjectType::SizeType radius )
    {
    numberOfTerms = 0;
    for( int i = 0; i < (int) images.size(); i++)
      {
      for( unsigned int j = 0; j < images[i].Size(); j++ )
        {
        if( images[i][j].GetPointer() )
          {
          faceLists.push_back( faceCalculator( images[i][j], region, radius ) );
          numberOfTerms++;
          }
        }
      }
    }

  void GoToBegin()
    {
    if( (int) faceListIts.size() != numberOfTerms )
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        faceListIts.push_back( faceLists[i].begin() );
        }
      }
    else
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        faceListIts[i] = faceLists[i].begin();
        }
      }
    }

  bool IsAtEnd()
    {
    for( int i = 0; i < numberOfTerms; i++ )
      {
      if( faceListIts[i] == faceLists[i].end() )
        {
        return true;
        }
      }
    return false;
    }

  void Increment()
    {
    for( int i = 0; i < numberOfTerms; i++ )
      {
      ++faceListIts[i];
      }
    }

  template< class IteratorType >
  void SetIteratorToCurrentFace(
      IteratorType& iterator,
      const ImageType& image,
      typename ImageType::ObjectType::SizeType radius )
    {
    if( image.GetPointer() )
      {
      iterator = IteratorType( radius, image, *faceListIts[0] );
      }
    else
      {
      iterator = IteratorType();
      }
    }

  template< class IteratorType >
  void SetIteratorToCurrentFace(
      IteratorType& iterator,
      const ImageType& image )
    {
    if( image.GetPointer() )
      {
      iterator = IteratorType( image, *faceListIts[0] );
      }
    else
      {
      iterator = IteratorType();
      }
    }

  template< class IteratorType >
  void SetIteratorToCurrentFace(
      std::vector< IteratorType >& iterators,
      const std::vector< ImageType >& images,
      typename ImageType::ObjectType::SizeType radius )
    {
    if( (int) iterators.size() != numberOfTerms )
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        if( images[i].GetPointer() )
          {
          iterators.push_back(
              IteratorType( radius, images[i], *faceListIts[i] ) );
          }
        else
          {
          iterators.push_back( IteratorType() );
          }
        }
      }
    else
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        if( images[i].GetPointer() )
          {
          iterators[i] = IteratorType( radius, images[i], *faceListIts[i] );
          }
        else
          {
          iterators[i] = IteratorType();
          }
        }
      }
    }

  template< class IteratorType >
  void SetIteratorToCurrentFace(
      std::vector< IteratorType >& iterators,
      const std::vector< ImageType >& images )
    {
    if( (int) iterators.size() != numberOfTerms )
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        if( images[i].GetPointer() )
          {
          iterators.push_back( IteratorType( images[i], *faceListIts[i] ) );
          }
        else
          {
          iterators.push_back( IteratorType() );
          }
        }
      }
    else
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        if( images[i].GetPointer() )
          {
          iterators[i] = IteratorType( images[i], *faceListIts[i] );
          }
        else
          {
          iterators[i] = IteratorType();
          }
        }
      }
    }

    template< class IteratorType, unsigned int VLength >
    void SetIteratorToCurrentFace(
        std::vector< itk::FixedArray< IteratorType, VLength > > &iterators,
        const std::vector< itk::FixedArray< ImageType, VLength > > & images )
      {
      int c = 0;
      if( (int) iterators.size() != (int) images.size() )
        {
        for( int i = 0; i < (int) images.size(); i++ )
          {
          itk::FixedArray< IteratorType, VLength > fixedArray;
          for( int j = 0; j < (int) images[i].Size(); j++ )
            {
            if( images[i][j] )
              {
              fixedArray[j] = IteratorType( images[i][j], *faceListIts[c] );
              c++;
              }
            else
              {
              fixedArray[j] = IteratorType();
              }
            }
          iterators.push_back( fixedArray );
          }
        }
      else
        {
        for( int i = 0; i < (int) images.size(); i++ )
          {
          itk::FixedArray< IteratorType, VLength > fixedArray;
          for( int j = 0; j < (int) images[i].Size(); j++ )
            {
            if( images[i][j].GetPointer() )
              {
              fixedArray[j] = IteratorType( images[i][j], *faceListIts[c] );
              c++;
              }
            else
              {
              fixedArray[j] = IteratorType();
              }
            }
          iterators[i] = fixedArray;
          }
        }
      }

  FaceCalculatorType                     faceCalculator;
  std::vector< FaceListType >            faceLists;
  std::vector< FaceListIteratorType >    faceListIts;
  int                                    numberOfTerms;
  };

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkDiffusiveRegistrationFilter+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkDiffusiveRegistrationFilter.txx"
#endif

#endif
