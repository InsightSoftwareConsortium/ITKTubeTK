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

#define itkDiffusiveRegistrationFilterNewMacro(x)              \
  static Pointer New(void)                                     \
    {                                                          \
    Pointer smartPtr = ::itk::ObjectFactory< x >::Create();    \
    if ( smartPtr.GetPointer() == NULL )                       \
      {                                                        \
      smartPtr = new x;                                        \
      }                                                        \
    smartPtr->UnRegister();                                    \
    smartPtr->CreateRegistrationFunction();                    \
    return smartPtr;                                           \
    }                                                          \
  virtual::itk::LightObject::Pointer CreateAnother(void) const \
    {                                                          \
    ::itk::LightObject::Pointer smartPtr;                      \
    smartPtr = x::New().GetPointer();                          \
    return smartPtr;                                           \
    }

namespace itk
{

/** \class itkDiffusiveRegistrationFilter
 * \brief Algorithm for registration of images depicting sliding organs, using
 * an anisotropic diffusive regularization term.
 *
 * Traditional deformable image registration imposes a uniform
 * smoothness constraint on the deformation field. This
 * is not appropriate when registering images visualizing organs
 * that slide relative to each other, and therefore leads to registration
 * inaccuracies.
 *
 * This algorithm includes a deformation field regularization term that is
 * based on anisotropic diffusion and accommodates the deformation field
 * discontinuities that are expected when considering sliding motion.
 *
 * The update term is composed of two parts: an intensity distance term and a
 * regularization term.  The intensity distance term uses the sum of square
 * difference metric, so this registration algorithm is appropriate for
 * monomodal image registration term only.  The regularization term uses a
 * specified border between the organs (stored as a vtkPolyData *) and enforces
 * coupling between the organs while allowing the motion field to exhibit
 * sliding motion at the organ interface.
 *
 * See: D.F. Pace et al., Deformable image registration of sliding organs using
 * anisotropic diffusive regularization, ISBI 2011.
 *
 * This class is templated over the type of the fixed image, the type of the
 * moving image and the type of the deformation field.
 *
 * \sa itkAnisotropicDiffusiveRegistrationFunction
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
  itkDiffusiveRegistrationFilterNewMacro(Self);

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
  typedef itk::ImageRegionIterator< OutputImageType > OutputImageRegionType;
  typedef typename FiniteDifferenceFunctionType::NeighborhoodType
      NeighborhoodType;
  typedef itk::ImageRegionIterator< UpdateBufferType > UpdateBufferRegionType;

  /** The registration function type */
  typedef DiffusiveRegistrationFunction
      < FixedImageType, MovingImageType, DeformationFieldType >
      RegistrationFunctionType;
  typedef typename RegistrationFunctionType::Pointer
      RegistrationFunctionPointer;
  typedef typename RegistrationFunctionType::RegularizationFunctionPointer
      RegularizationFunctionPointer;
  typedef typename RegistrationFunctionType::SpacingType    SpacingType;

  /** Deformation field types. */
  typedef typename RegistrationFunctionType::DeformationVectorType
      DeformationVectorType;
  typedef typename RegistrationFunctionType::DeformationVectorComponentType
      DeformationVectorComponentType;
  typedef typename RegistrationFunctionType::DeformationVectorComponentImageType
      DeformationVectorComponentImageType;
  typedef typename DeformationVectorComponentImageType::Pointer
      DeformationVectorComponentImagePointer;
  typedef typename
      RegistrationFunctionType::DeformationVectorComponentNeighborhoodType
      DeformationVectorComponentNeighborhoodType;
  typedef typename
      RegistrationFunctionType::DeformationVectorComponentNeighborhoodArrayType
      DeformationVectorComponentNeighborhoodArrayType;
  typedef typename DeformationVectorComponentImageType::RegionType
      ThreadDeformationVectorComponentImageRegionType;

  /** Diffusion tensor image types */
  typedef typename RegistrationFunctionType::DiffusionTensorImageType
      DiffusionTensorImageType;
  typedef typename DiffusionTensorImageType::Pointer
      DiffusionTensorImagePointer;
  typedef typename RegistrationFunctionType::DiffusionTensorNeighborhoodType
      DiffusionTensorNeighborhoodType;
  typedef typename DiffusionTensorImageType::RegionType
      ThreadDiffusionTensorImageRegionType;

  /** Tensor derivative matrix image types */
  typedef typename RegistrationFunctionType::TensorDerivativeImageType
      TensorDerivativeImageType;
  typedef typename TensorDerivativeImageType::Pointer
      TensorDerivativeImagePointer;
  typedef typename RegistrationFunctionType::TensorDerivativeImageRegionType
      TensorDerivativeImageRegionType;
  typedef typename TensorDerivativeImageType::RegionType
      ThreadTensorDerivativeImageRegionType;

  /** Deformation component image array types */
  typedef typename
      itk::FixedArray< DeformationVectorComponentImagePointer, ImageDimension >
      DeformationComponentImageArrayType;

  /** Types for vector component extractor */
  typedef itk::VectorIndexSelectionCastImageFilter
      < DeformationFieldType, DeformationVectorComponentImageType >
      VectorIndexSelectionFilterType;

  /** Convenience functions to set/get the registration functions timestep. */
  void SetTimeStep( const TimeStepType &t )
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

  /** Get the image of the diffusion tensor */
  virtual const DiffusionTensorImagePointer GetDiffusionTensorImage() const
    { return m_DiffusionTensorImage; }

  /** Get the image of the diffusion tensor derivatives */
  virtual const TensorDerivativeImagePointer GetDiffusionTensorDerivativeImage()
      const
    { return m_DiffusionTensorDerivativeImage; }

protected:
  DiffusiveRegistrationFilter();
  virtual ~DiffusiveRegistrationFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Initialization occuring before the registration iterations. */
  virtual void Initialize();

  /** Allocate images used during the registration. */
  virtual void AllocateImages();

  /** Computes the diffusion tensor images */
  virtual void ComputeDiffusionTensorImages();

  /** Computes the first derivatives of the diffusion tensor images */
  virtual void ComputeDiffusionTensorDerivativeImages();

  /** Helper to compute the first derivatives of the diffusion tensor images */
  virtual void ComputeDiffusionTensorDerivativeImageHelper(
      DiffusionTensorImagePointer tensorImage,
      TensorDerivativeImagePointer tensorDerivativeImage );

  /** Allocate the update buffer. */
  virtual void AllocateUpdateBuffer();

  /** Get the update buffer. */
  virtual UpdateBufferType * GetUpdateBuffer()
    { return m_UpdateBuffer; }

  /** Initialize the state of the filter and equation before each iteration. */
  virtual void InitializeIteration();

  /** Updates the deformation vector component images */
  virtual void UpdateDeformationVectorComponentImages();

  /** Extracts the x, y, z components of a deformation field. */
  void ExtractXYZComponentsFromDeformationField(
      OutputImagePointer deformationField,
      DeformationComponentImageArrayType& deformationComponentImages );

  /** Get the array of deformation component images. */
  DeformationComponentImageArrayType GetDeformationComponentImages()
    { return m_DeformationComponentImages; }

  /** This method populates an update buffer with changes for each pixel in the
   * output, using the ThreadedCalculateChange() method and a multithreading
   * mechanism. Return value is a time step to be used for the update.
   * \sa ThreadedCalculateChange */
  virtual TimeStepType CalculateChange();

  /** Inherited from superclass - do not call this function!  Call the other
   *  ThreadedCalculateChange function instead */
  TimeStepType ThreadedCalculateChange( const ThreadRegionType &regionToProcess,
                                       int threadId );

  /** Does the actual work of calculating change over a region supplied by
   * the multithreading mechanism.
   * \sa CalculateChange
   * \sa CalculateChangeThreaderCallback */
  virtual TimeStepType ThreadedCalculateChange(
      const ThreadRegionType &regionToProcess,
      const ThreadDiffusionTensorImageRegionType &tensorRegionToProcess,
      const ThreadTensorDerivativeImageRegionType
        &tensorDerivativeRegionToProcess,
      const ThreadDeformationVectorComponentImageRegionType
        &deformationComponentRegionToProcess,
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
                                    const ThreadRegionType &regionToProcess,
                                    int threadId );

  /** Helper function to allocate an image based on a template */
  template< class UnallocatedImageType, class TemplateImageType >
  void AllocateSpaceForImage( UnallocatedImageType & image,
                             const TemplateImageType & templateImage );

  /** Helper function to check whether the attributes of an image match a
    * template */
  template< class CheckedImageType, class TemplateImageType >
  bool CompareImageAttributes( const CheckedImageType & image,
                               const TemplateImageType & templateImage );

  /** Create the registration function, with default parameters for
    * ComputeRegularizationTerm and ComputeIntensityDistanceTerm. */
  virtual void CreateRegistrationFunction();

  /** Get the registration function pointer */
  virtual RegistrationFunctionType * GetRegistrationFunctionPointer() const;

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

  /** Image storing information we will need for each voxel on every
   *  registration iteration */
  DiffusionTensorImagePointer         m_DiffusionTensorImage;
  TensorDerivativeImagePointer        m_DiffusionTensorDerivativeImage;
  DeformationComponentImageArrayType  m_DeformationComponentImages;
};

/** Struct to simply get the face list and an iterator over the face list
 *  when processing an image */
template< class ImageType >
struct FaceStruct
  {
  FaceStruct() {}

  FaceStruct( ImageType& image,
              typename ImageType::ObjectType::SizeType radius )
    {
    if( image )
      {
      faceList = faceCalculator( image,
                                 image->GetLargestPossibleRegion(),
                                 radius );
      }
    }

  FaceStruct( ImageType& image,
              typename ImageType::ObjectType::RegionType region,
              typename ImageType::ObjectType::SizeType radius )
    {
    if( image )
      {
      faceList = faceCalculator( image, region, radius );
      }
    }

  void begin()
    {
    faceListIt = faceList.begin();
    }

  bool IsAtEnd()
    {
    return faceListIt == faceList.end();
    }

  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
      < typename ImageType::ObjectType > FaceCalculatorType;
  typedef typename FaceCalculatorType::FaceListType FaceListType;
  typedef typename FaceListType::iterator FaceListIteratorType;

  FaceCalculatorType      faceCalculator;
  FaceListType            faceList;
  FaceListIteratorType    faceListIt;
  };

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkDiffusiveRegistrationFilter+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkDiffusiveRegistrationFilter.txx"
#endif

#endif
