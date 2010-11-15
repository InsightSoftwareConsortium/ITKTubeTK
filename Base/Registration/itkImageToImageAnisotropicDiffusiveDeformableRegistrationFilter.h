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
#ifndef __itkImageToImageAnisotropicDiffusiveDeformableRegistrationFilter_h
#define __itkImageToImageAnisotropicDiffusiveDeformableRegistrationFilter_h

#include "itkPDEDeformableRegistrationFilter.h"
#include "itkImageToImageAnisotropicDiffusiveDeformableRegistrationFunction.h"

#include "vtkPointLocator.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"

namespace itk
{

/** \class itkImageToImageAnisotropicDiffusiveDeformableRegistrationFilter
 * \brief Insert brief description here!!!
 *
 * Implements diffusive deformable registration, where the update term
 * is composed of two parts: an intensity difference term and a regularization
 * term.  The intensity difference term is computed based on sum of square
 * differences, with the assumption that this algorithm will be used for
 * monomodal image registration at the moment.  The regularization term ensures
 * that the resulting deformation field is realistic, and has been tailored for
 * the registration of images depicting images that slide relative to each
 * other.  By specifying the border between the organs (using a vtkPolyData * )
 * we can ensure that the motion field is smooth in the direction parallel to
 * the border's normal (to enforce coupling between the organs) but allow
 * the motion field to be discontinuous in the direction parallel to the border
 * itself in the vicinity of the border (to allow for sliding motion).
 *
 * Insert paper reference here!!!!!!!
 *
 * Insert more description + warnings here!!!!!
 *
 * This class is templated over the fixed image type, moving image type and the
 * deformation field type.
 *
 * \sa itkImageToImageAnisotropicDiffusiveDeformableRegistrationFunction
 * \ingroup DeformableImageRegistration MultiThreaded
 */

template < class TFixedImage, class TMovingImage, class TDeformationField >
class ITK_EXPORT ImageToImageAnisotropicDiffusiveDeformableRegistrationFilter
  : public PDEDeformableRegistrationFilter< TFixedImage,
                                           TMovingImage,
                                           TDeformationField >
{
public:
  /** Standard class typedefs. */
  typedef ImageToImageAnisotropicDiffusiveDeformableRegistrationFilter
        Self;
  typedef PDEDeformableRegistrationFilter< TFixedImage,
                                          TMovingImage,
                                          TDeformationField > Superclass;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, PDEDeformableRegistrationFilter);

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Inherit FiniteDifferenceFunction type from the superclass. */
  typedef typename Superclass::FiniteDifferenceFunctionType
      FiniteDifferenceFunctionType;

  /** Types for the fixed image. */
  typedef typename Superclass::FixedImageType           FixedImageType;
  typedef typename Superclass::FixedImagePointer        FixedImagePointer;

  /** Types for the moving image. */
  typedef typename Superclass::MovingImageType          MovingImageType;
  typedef typename Superclass::MovingImagePointer       MovingImagePointer;

  /** Types for the deformation field. */
  typedef typename Superclass::DeformationFieldType     DeformationFieldType;
  typedef typename Superclass::DeformationFieldPointer  DeformationFieldPointer;

  /** Inherit types from superclass. */
  typedef typename Superclass::TimeStepType             TimeStepType;

  /** The registration function type */
  typedef ImageToImageAnisotropicDiffusiveDeformableRegistrationFunction
      < FixedImageType, MovingImageType, DeformationFieldType >
      RegistrationFunctionType;

  /** Deformation field types - types for the deformation vectors, deformation
   *  vector components, and vector component images
   */
  typedef typename RegistrationFunctionType::DeformationVectorType
      DeformationVectorType;
  typedef typename RegistrationFunctionType::DeformationVectorComponentType
      DeformationVectorComponentType;
  typedef typename RegistrationFunctionType::DeformationVectorComponentImageType
      DeformationVectorComponentImageType;
  typedef typename DeformationVectorComponentImageType::Pointer
      DeformationVectorComponentImagePointer;

  /** Normal vector types - types for the normal vectors and normal vector
   *  images */
  typedef typename RegistrationFunctionType::NormalVectorType
      NormalVectorType;
  typedef typename RegistrationFunctionType::NormalVectorImageType
      NormalVectorImageType;
  typedef typename NormalVectorImageType::Pointer
      NormalVectorImagePointer;
  typedef itk::ImageRegionIterator< NormalVectorImageType >
      NormalVectorImageIteratorType;

  /** Weight image types - types for the weightings and weight images.  This
    * is the weight between the anisotropic diffusive regularization and a
    * diffusive (Gaussian) regularization.
    */
  typedef double                                        WeightType;
  typedef itk::Image< WeightType, ImageDimension >      WeightImageType;
  typedef typename WeightImageType::Pointer             WeightImagePointer;
  typedef itk::ImageRegionIterator< WeightImageType >   WeightImageIteratorType;

  /** Organ boundary surface types - types for the surface and the surface of
    * normals, and point locator to find points on the surface */
  typedef vtkSmartPointer< vtkPolyData >                BorderSurfacePointer;
  typedef vtkPolyDataNormals
      BorderNormalsSurfaceFilterType;
  typedef vtkSmartPointer< BorderNormalsSurfaceFilterType >
      BorderNormalsSurfaceFilterPointer;
  typedef vtkPointLocator                               PointLocatorType;
  typedef vtkSmartPointer< PointLocatorType >           PointLocatorPointer;

  /** The diffusion tensor types */
  typedef typename RegistrationFunctionType::DiffusionTensorImageType
      DiffusionTensorImageType;
  typedef typename RegistrationFunctionType::DiffusionTensorImagePointer
      DiffusionTensorImagePointer;

  /** Typedefs used in multithreading */
  typedef typename Superclass::OutputImageType          OutputImageType;
  typedef typename Superclass::OutputImagePointer       OutputImagePointer;
  typedef typename Superclass::UpdateBufferType         UpdateBufferType;

  /** Region types used in multithreading */
  typedef typename UpdateBufferType::RegionType         ThreadRegionType;
  typedef typename NormalVectorImageType::RegionType
      ThreadNormalVectorImageRegionType;
  typedef typename DiffusionTensorImageType::RegionType
      ThreadDiffusionTensorImageRegionType;
  typedef typename DeformationVectorComponentImageType::RegionType
      ThreadDeformationVectorComponentImageRegionType;

  /** Types for vector component extractor */
  typedef itk::VectorIndexSelectionCastImageFilter
      < DeformationFieldType, DeformationVectorComponentImageType >
      SelectionCastImageFilterType;
  typedef typename SelectionCastImageFilterType::Pointer
      SelectionCastImageFilterPointer;

  /** Convenience function to set the registration function's timestep from the
   * filter */
  void SetTimeStep( const TimeStepType &t );
  const TimeStepType& GetTimeStep() const;

  /** Set/get the organ boundary polydata, which must be in the same space as
   *  the fixed image.  Border normals are computed based on this polydata/
   */
  virtual void SetBorderSurface( BorderSurfacePointer border )
    { m_BorderSurface = border; }
  virtual const BorderSurfacePointer GetBorderSurface() const
    { return m_BorderSurface; }

  /** Get the polydata containing the border surface normals */
  virtual const BorderSurfacePointer GetBorderNormalsSurface() const
    { return m_BorderNormalsSurface; }

  /** Set/get the lambda that controls the exponential decay used to calculate
   *  the weight value w as a function of the distance to the closest border
   *  point.  Must be negative. */
  void SetLambda( WeightType l )
    { if ( l < 0 ) { m_lambda = l; } }
  WeightType GetLambda() const
    { return m_lambda; }

  /** Set/get whether to compute the intensity distance term
   *  Default: true
   */
  void SetComputeIntensityDistanceTerm( bool compute );
  bool GetComputeIntensityDistanceTerm() const;

  /** Set/get whether to compute the motion field regularization term
   *  Default: true
   */
  void SetComputeRegularizationTerm( bool compute );
  bool GetComputeRegularizationTerm() const;

  /** Set/get whether to use the anisotropic diffusive regularization.  If
   *  false, the weighting term w=0 and Gaussian regularization is used.
   *  Default: true
   */
  void SetUseAnisotropicRegularization( bool diffuse )
    { m_UseAnisotropicRegularization = diffuse; }
  bool GetUseAnisotropicRegularization() const
    { return m_UseAnisotropicRegularization; }

  /** Get/get the image of the normal vectors.  Setting the normal vector image
    * overrides the border surface polydata if a border surface was also
    * supplied. */
  virtual void SetNormalVectorImage( NormalVectorImagePointer normalImage )
    { m_NormalVectorImage = normalImage; }
  virtual const NormalVectorImagePointer GetNormalVectorImage() const
    { return m_NormalVectorImage; }

  /** Set/get the weighting image.  Setting the weighting image overrides
    * the border surface polydata if a border surface was also supplied.
    */
  virtual void SetWeightImage( WeightImagePointer weightImage )
    { m_WeightImage = weightImage; }
  virtual const WeightImagePointer GetWeightImage() const
    { return m_WeightImage; }

  /** Get the image of the tangential diffusion tensors */
  virtual const DiffusionTensorImagePointer GetTangentialDiffusionTensorImage()
      const
    { return m_TangentialDiffusionTensorImage; }

  /** Get the image of the normal diffusion tensors */
  virtual const DiffusionTensorImagePointer GetNormalDiffusionTensorImage()
      const
    { return m_NormalDiffusionTensorImage; }

  /** Get the normal components of the deformation field */
  virtual const OutputImagePointer GetNormalDeformationFieldImage() const
    { return m_NormalDeformationField; }

protected:
  ImageToImageAnisotropicDiffusiveDeformableRegistrationFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Initialization occuring before the registration loop. */
  virtual void Initialize();

  /** Allocate the update buffer - reimplented here to also allocate other
   *  images used in the registration computation. */
  virtual void AllocateUpdateBuffer();

  /** Initialize the state of the filter and equation before each iteration. */
  virtual void InitializeIteration();

  /** This method populates an update buffer with changes for each pixel in the
   * output, using the ThreadedCalculateChange() method and a multithreading
   * mechanism. Return value is a time step to be used for the update.
   * \sa ThreadedCalculateChange */
  virtual TimeStepType CalculateChange();

  /** Inherited from superclass - do not call this function!  Call the other
   *  ThreadedCalculateChange instead */
  TimeStepType ThreadedCalculateChange(
      const ThreadRegionType & regionToProcess, int threadId );

  /** Does the actual work of calculating change over a region supplied by
   * the multithreading mechanism.
   * \sa CalculateChange
   * \sa CalculateChangeThreaderCallback */
  virtual
  TimeStepType ThreadedCalculateChange(
          const ThreadRegionType &regionToProcess,
          const ThreadNormalVectorImageRegionType
            &normalVectorRegionToProcess,
          const ThreadDiffusionTensorImageRegionType
            &diffusionRegionToProcess,
          const ThreadDeformationVectorComponentImageRegionType
            &componentRegionToProcess,
          int threadId );

  /** This method applies changes from the update buffer to the output, using
   * the ThreadedApplyUpdate() method and a multithreading mechanism.  "dt" is
   * the time step to use for the update of each pixel.
   * \sa ThreadedApplyUpdate */
  virtual void ApplyUpdate(TimeStepType dt);

  /**  Does the actual work of updating the output from the UpdateContainer over
   *  an output region supplied by the multithreading mechanism.
   *  \sa ApplyUpdate
   *  \sa ApplyUpdateThreaderCallback */
  virtual
  void ThreadedApplyUpdate(
      TimeStepType dt, const ThreadRegionType &regionToProcess, int threadId );

  /** Computes the normal vector image and weighting factors w given the
   *  surface border polydata.
   */
  virtual void ComputeNormalVectorAndWeightImages(
      bool computeNormalVectorImage, bool computeWeightImage );

  /** Computes the weighting factor w from the distance to the border.  The
   *  weight should be 1 near the border and 0 away from the border.
   */
  virtual WeightType ComputeWeightFromDistance( WeightType distance );

  /** Updates the deformation vector component images */
  virtual void UpdateDeformationVectorComponentImages();

  /** Computes the diffusion tensor image */
  virtual void ComputeDiffusionTensorImage();

  /** Helper function to allocate an image based on a template */
  template< class UnallocatedImageType, class TemplateImageType >
  void AllocateSpaceForImage(
      UnallocatedImageType& inputImage, const TemplateImageType& templateImage );

  /** Helper function to check whether the attributes of an image match a
    * template */
  template< class CheckedImageType, class TemplateImageType >
  bool CompareImageAttributes(
      const CheckedImageType& inputImage, const TemplateImageType& templateImage);

private:
  // Purposely not implemented
  ImageToImageAnisotropicDiffusiveDeformableRegistrationFilter(const Self&);
  void operator=(const Self&); // Purposely not implemented

  /** Structure for passing information into static callback methods.  Used in
   * the subclasses' threading mechanisms. */
  struct DenseFDThreadStruct
    {
    ImageToImageAnisotropicDiffusiveDeformableRegistrationFilter *Filter;
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
  typename UpdateBufferType::Pointer    m_UpdateBuffer;

  /** Whether or not to use the diffusive regularization. */
  bool                                  m_UseAnisotropicRegularization;

  /** The organ boundary surface, the surface of border normals, and the derived
   *  normal vector and weight images */
  BorderSurfacePointer                  m_BorderSurface;
  BorderSurfacePointer                  m_BorderNormalsSurface;
  BorderNormalsSurfaceFilterPointer     m_BorderNormalsSurfaceFilter;
  NormalVectorImagePointer              m_NormalVectorImage;
  WeightImagePointer                    m_WeightImage;

  /** Computes the nearest point on the polydata to a given point */
  PointLocatorPointer                   m_PointLocator;

  /** The lambda factor for computing the weight from distance.  Weight is
    * modeled as exponential decay: weight = e^(lambda * distance).
    * (lamba must be negative)
    */
  WeightType                            m_lambda;

  /** The normal component of the deformation field */
  OutputImagePointer                    m_NormalDeformationField;

  /** The components of the tangential and normal deformation vectors */
  itk::FixedArray< DeformationVectorComponentImagePointer, ImageDimension >
      m_DeformationVectorTangentialComponents;
  itk::FixedArray< DeformationVectorComponentImagePointer, ImageDimension >
      m_DeformationVectorNormalComponents;

  /** Extracts the x,y,z components of the tangential and normal components of
   * the deformation field
   */
  itk::FixedArray< SelectionCastImageFilterPointer, ImageDimension >
                                        m_TangentialComponentExtractor;
  itk::FixedArray< SelectionCastImageFilterPointer, ImageDimension >
                                        m_NormalComponentExtractor;

  /** The images of the tangential and normal diffusion tensors */
  DiffusionTensorImagePointer           m_TangentialDiffusionTensorImage;
  DiffusionTensorImagePointer           m_NormalDiffusionTensorImage;

};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkImageToImageAnisotropicDiffusiveDeformableRegistrationFilter+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkImageToImageAnisotropicDiffusiveDeformableRegistrationFilter.txx"
#endif

#endif

/** TODO LIST - essential

  // Calculate w in ComputeNormalVectorAndWeightImages() - but don't bother if
  // m_UseAnisotropicRegularization == false (do it in the loop with the normals)

  // Get w in ComputeDiffusionTensorImage()

  // checking the timestep for stability as in the anisotropic filter

  // halting criteria?!?!

  // SetUseImageSpacing() to on for this filter?  Would give derivates in physical space, default is off

  // why did Andinet comment out m_UpdateBuffer->Modified() in CalculateChange()

  // boundary faces

  // difference between cell normals and normals - using the correct ones?

  // better way to compute n's?

  // Better (faster?) alternative to point locator

  // Go through PDERegistrationFilter and PDERegistrationFunction classes to
  // ensure we're not missing any important functions to override / parameters
  // to set

  // see where you can make it faster - computing anything unnecessarily?

  ======== function =========

  // test to see if andinet's stuff works with a timestep of 1.0

  // compute timestep instead of hard-coding 0.05 in constructor

  // does the intensiy distance function need more parameters in the ctor?
  // Ex. options for the global filter
  // ... inherited for finite
  //  m_IntensityDistanceFunction->SetScaleCoefficients( vals );
  // ... inherited from PDE function
  //  m_IntensityDistanceFunction->SetGradientStep( 0.0 );
  //  m_IntensityDistanceFunction->SetNormalizeGradient( false );

  // In ctor - intensity distance function uses LinearInterpolateImageFilter
  // by default, can change later using setMovingImageInterpolator()

  // check to make sure it's ok that intensity function uses global timestep
  // i think it is because we are just calling update(), which does not use
  // m_TimeStep

  // Better std::cout for what's being done -with iteration number

  // make sure that there is no normalization in the intensity distance
  // function - check that smooth gradient is off, and that the normalize
  // metric doesn't do anything we don't want

    // weighting between the intensity distance and regularization terms
  // in ComputeUpdate?  Don't worry about weighting if one term is not
  // computed because of boolean settings


  */

/** TODO LIST - medium

  // Allocating all of the images ( in Initialize() ) essential?

    ======== function =========

  */

/** TODO LIST - later

  // Better class description at top

  // go through typedefs and pull as many as possible from the function!
  // may be able to simplify the typedefs here to make them clearer

  // Specification of vtkPolyData introduces dependency on VTK and is
  // horrible for dimension-independence... perhaps provide additional /
  // only option to provide normal image only - then could ensure that it
  // with the right templates

  // printSelf not compiling because of os type mismatch

  // create superclass with these methods for anisotropic diffusion
  // registration filter

  // assert vs. if?

  // progress object as in demons test

  // pixel type vs. deformation field type?  deformation field type should
  // always be double?

  // what happens if normals are not set to registrator? i.e. normals are
  // the defaults of 0,0,0

    ======== function =========

    // Better class description at top

    // ComputeUpdate() won't work, but shouldn't be called anyways

    ====== tests =======

    // try it out with dimension = 2

  // there are some exception handling tests in
  // itkDemonsRegistrationFilterTest that would be good to put in the image
    // registration test

   // experiment with adding noise to the border of the motion field
    // regularization test

  */
