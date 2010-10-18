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
#ifndef __itkImageToImageDiffusiveDeformableRegistrationFilter_h
#define __itkImageToImageDiffusiveDeformableRegistrationFilter_h

#include "itkPDEDeformableRegistrationFilter.h"
#include "itkImageToImageDiffusiveDeformableRegistrationFunction.h"

namespace itk
{

/** \class itkImageToImageDiffusiveDeformableRegistrationFilter
 * \brief TODO
 *
 * TODO Insert more description + warnings here
 *
 * This class is templated over the fixed image type, moving image type and the
 * deformation field type.
 *
 * \sa itkImageToImageDiffusiveDeformableRegistrationFunction
 * \ingroup DeformableImageRegistration MultiThreaded
 */

template < class TFixedImage, class TMovingImage, class TDeformationField >
class ITK_EXPORT ImageToImageDiffusiveDeformableRegistrationFilter
  : public PDEDeformableRegistrationFilter< TFixedImage,
                                           TMovingImage,
                                           TDeformationField >
{
public:
  /** Standard class typedefs. */
  typedef ImageToImageDiffusiveDeformableRegistrationFilter   Self;
  typedef PDEDeformableRegistrationFilter< TFixedImage,
                                          TMovingImage,
                                          TDeformationField > Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, PDEDeformableRegistrationFilter);

  /** Inherit types from superclass. */
  typedef typename Superclass::TimeStepType             TimeStepType;

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::FixedImagePointer        FixedImagePointer;

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType          MovingImageType;
  typedef typename Superclass::MovingImagePointer       MovingImagePointer;

  /** Deformation field type. */
  typedef typename Superclass::DeformationFieldType     DeformationFieldType;
  typedef typename Superclass::DeformationFieldPointer  DeformationFieldPointer;

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  // TODO go through typedefs and pull as many as possible from the function!
  // TODO may be able to simplify the typedefs here to make them clearer

  /** Deformation field types. */

  // ex. vector < double, 3 >
  typedef typename Superclass::DeformationFieldType::PixelType
                                    DeformationFieldVectorType;
  // ex. double
  typedef typename DeformationFieldType::PixelType::ValueType
                                    DeformationFieldScalarType;
  // ex. image of doubles
  typedef itk::Image< DeformationFieldScalarType, ImageDimension >
                                    DeformationFieldComponentImageType;
  typedef typename DeformationFieldComponentImageType::Pointer
                                    DeformationFieldComponentImagePointer;
  typedef itk::ImageRegionIterator< DeformationFieldComponentImageType >
                                    DeformationFieldComponentImageIteratorType;

  /** Normal vector types. */
  typedef DeformationFieldVectorType                     NormalVectorType;
  typedef itk::Image< NormalVectorType, ImageDimension > NormalVectorImageType;
  typedef typename NormalVectorImageType::Pointer        NormalVectorImagePointer;
  typedef itk::ImageRegionIterator< NormalVectorImageType >
                                                         NormalVectorIteratorType;

  /** FiniteDifferenceFunction type. */
  typedef typename Superclass::FiniteDifferenceFunctionType
                                                  FiniteDifferenceFunctionType;

  /** ImageToImageDiffusiveDeformableRegistrationFunction type. */
  typedef ImageToImageDiffusiveDeformableRegistrationFunction
                                            < FixedImageType,
                                              MovingImageType,
                                              DeformationFieldType >
                                              RegistrationFunctionType;
  typedef typename RegistrationFunctionType::DiffusionTensorImageType
                                              DiffusionTensorImageType;
  typedef typename DiffusionTensorImageType::PixelType
                                              DiffusionTensorImagePixelType;
  typedef typename DiffusionTensorImageType::Pointer
                                              DiffusionTensorImagePointer;
  typedef typename RegistrationFunctionType::DefaultBoundaryConditionType
                                              DefaultBoundaryConditionType;

  /** The type of region used for multithreading */
  typedef typename Superclass::OutputImageType  OutputImageType;
  typedef OutputImageType                       UpdateBufferType;
  typedef typename UpdateBufferType::RegionType ThreadRegionType;

  /** The type of region used for multithreading */
  typedef typename NormalVectorImageType::RegionType
                                ThreadNormalVectorImageRegionType;
  typedef typename DiffusionTensorImageType::RegionType
                                ThreadDiffusionTensorImageRegionType;
  typedef typename DeformationFieldComponentImageType::RegionType
                                ThreadDeformationFieldComponentImageRegionType;

  /** Define neighborhood types */
  typedef typename RegistrationFunctionType::NormalVectorImageNeighborhoodType
                                NormalVectorImageNeighborhoodType;

  typedef ConstNeighborhoodIterator< DiffusionTensorImageType,
                                DefaultBoundaryConditionType >
                                DiffusionTensorNeighborhoodType;

  typedef typename
        RegistrationFunctionType::DeformationFieldComponentNeighborhoodType
                                DeformationFieldComponentNeighborhoodType;
  typedef typename
        RegistrationFunctionType::DeformationFieldComponentNeighborhoodArrayType
                                DeformationFieldComponentNeighborhoodArrayType;

  /** Set the border normal. The magnitude of the vector is internally set
   *  to one if it isn't already.
   */
  // TODO need to calculate this here or input as image of normals.  For now
  // we're taking only one vector for this test.
  virtual void SetNormalVectors( NormalVectorType& normals );
  virtual const NormalVectorType& GetNormalVectors() const;

  /** Get the image of the tangential diffusion tensors */
  virtual const DiffusionTensorImagePointer GetTangentialDiffusionTensorImage()
                                                                        const;

  /** Get the image of the normal diffusion tensors */
  virtual const DiffusionTensorImagePointer GetNormalDiffusionTensorImage()
                                                                        const;

  /** Whether to compute the motion field regularization term (for testing)
   *  Default: true */
  void SetComputeRegularizationTerm( bool compute );
  bool GetComputeRegularizationTerm() const;

  /** Whether to compute the intensity distance term (for testing)
   *  Default: true */
  void SetComputeIntensityDistanceTerm( bool compute );
  bool GetComputeIntensityDistanceTerm() const;

protected:
  ImageToImageDiffusiveDeformableRegistrationFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Initialize the state of the filter and equation before each iteration. */
  virtual void InitializeIteration();

  /** This method populates an update buffer with changes for each pixel in the
   * output using the ThreadedCalculateChange() method and a multithreading
   * mechanism. Returns value is a time step to be used for the update. */
  virtual TimeStepType CalculateChange();

  /** Allocate the deformation field component images */
  virtual void AllocateDeformationFieldComponentImages();

  /** Allocate the diffusion tensor image. */
  virtual void AllocateDiffusionTensorImage();

  /** Allocate the update buffer - reimplemented here to also call
   *  AllocateDiffusionTensorImage. */
  virtual void AllocateUpdateBuffer();

  /** Helper function to allocate an image based on a template */
  template< class UnallocatedImageType, class TemplateImageType >
  void AllocateSpaceForImage( UnallocatedImageType& inputImage,
                              const TemplateImageType& templateImage );

  /** Update the normal vector image and weighting factor w */
  virtual void UpdateNormalVectorImage();

  /** Update diffusion tensor image */
  virtual void UpdateDiffusionTensorImage();

  /** Update deformation field component images */
  virtual void UpdateDeformationFieldComponentImages();

  /** This method applies changes from the m_UpdateBuffer to the output using
   * the ThreadedApplyUpdate() method and a multithreading mechanism.  "dt" is
   * the time step to use for the update of each pixel. */
  virtual void ApplyUpdate(TimeStepType dt);

  /**  Does the actual work of updating the output from the UpdateContainer over
   *  an output region supplied by the multithreading mechanism.
   *  \sa ApplyUpdate
   *  \sa ApplyUpdateThreaderCallback */
  virtual
  void ThreadedApplyUpdate(TimeStepType dt,
                           const ThreadRegionType &regionToProcess,
                           const ThreadDiffusionTensorImageRegionType
                                                  &diffusionTensorImageRegion,
                           int threadId);

  /** Does the actual work of calculating change over a region supplied by
   * the multithreading mechanism.
   * \sa CalculateChange
   * \sa CalculateChangeThreaderCallback */
  virtual
  TimeStepType ThreadedCalculateChange(
          const ThreadRegionType &regionToProcess,
          const ThreadNormalVectorImageRegionType &normalVectorRegionToProcess,
          const ThreadDiffusionTensorImageRegionType &diffusionRegionToProcess,
          const ThreadDeformationFieldComponentImageRegionType
                                                      &componentRegionToProcess,
          int threadId);

private:
  // Purposely not implemented
  ImageToImageDiffusiveDeformableRegistrationFilter(const Self&);
  void operator=(const Self&); // Purposely not implemented

  /** Structure for passing information into static callback methods.  Used in
   * the subclasses' threading mechanisms. */
  struct DenseFDThreadStruct
    {
    ImageToImageDiffusiveDeformableRegistrationFilter *Filter;
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

  /** The buffer that holds the updates for an iteration of the algorithm. */
  typename UpdateBufferType::Pointer m_UpdateBuffer;

  /** The border normals. */
  // TODO take me out
  NormalVectorType                      m_NormalVectors;
  NormalVectorImagePointer              m_NormalVectorImage;

  /** The image of the tangential diffusion tensors
   * Calculate div( P^TP gradient(motion field tangential component image) )
   */
  DiffusionTensorImagePointer           m_TangentialDiffusionTensorImage;

  /** The image of the normal diffusion tensors
   * Calculate div( w(n^T gradient(motion field normal component image) ) n )
   */
  DiffusionTensorImagePointer           m_NormalDiffusionTensorImage;

  /** Extracts the tangential and normal components of the deformation field */
  typedef typename OutputImageType::Pointer OutputImagePointer;
  OutputImagePointer                    m_OutputTangentialImage;
  OutputImagePointer                    m_OutputNormalImage;

  /** Extracts the x,y,z components of the tangential and normal components of
   * the deformation field
   */
  typedef itk::VectorIndexSelectionCastImageFilter< DeformationFieldType,
                                        DeformationFieldComponentImageType >
                                        SelectionCastImageFilterType;
  typedef typename SelectionCastImageFilterType::Pointer
                                        SelectionCastImageFilterPointer;
  itk::FixedArray< SelectionCastImageFilterPointer, ImageDimension >
                                        m_TangentialComponentExtractor;
  itk::FixedArray< SelectionCastImageFilterPointer, ImageDimension >
                                        m_NormalComponentExtractor;

  itk::FixedArray< DeformationFieldComponentImagePointer, ImageDimension >
                                        m_DeformationFieldTangentialComponents;

  itk::FixedArray< DeformationFieldComponentImagePointer, ImageDimension >
                                        m_DeformationFieldNormalComponents;

};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkImageToImageDiffusiveDeformableRegistrationFilter+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkImageToImageDiffusiveDeformableRegistrationFilter.txx"
#endif

#endif
