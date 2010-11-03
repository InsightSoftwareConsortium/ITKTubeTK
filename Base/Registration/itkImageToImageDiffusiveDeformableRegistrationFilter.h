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

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataNormals.h"
#include "vtkPointLocator.h"

namespace itk
{

/** \class itkImageToImageDiffusiveDeformableRegistrationFilter
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

  /** Weight image types. */
  typedef DeformationFieldScalarType                    WeightType;
  typedef itk::Image< WeightType, ImageDimension >      WeightImageType;
  typedef typename WeightImageType::Pointer             WeightImagePointer;
  typedef itk::ImageRegionIterator< WeightImageType >   WeightIteratorType;

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

  /** Convenience function to set the function's timestep from the filter */
  void SetTimeStep( const TimeStepType &t );
  const TimeStepType& GetTimeStep() const;

  /** The type of region used for multithreading */
  typedef typename Superclass::OutputImageType  OutputImageType;
  typedef typename OutputImageType::Pointer OutputImagePointer;
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

  /** Set the border polydata for the fixed image.  Border normals are computed
   *  based on this polydata, so it should be "well-behaved" under
   *  vtkPolyDataNormals.  The border normal must be in the same space as the
   *  fixed image.
   */
  typedef vtkPolyData                           BorderSurfaceType;
  typedef vtkSmartPointer< BorderSurfaceType >  BorderSurfacePointer;
  virtual void SetBorderSurface( BorderSurfacePointer border );
  virtual const BorderSurfacePointer GetBorderSurface() const
    { return m_BorderSurface; }

  /** Get the polydata of the border surface normals */
  virtual const BorderSurfacePointer GetBorderNormalsSurface() const
    { return m_BorderNormalsSurface; }

  /** Get the image of the normal vectors */
  virtual const NormalVectorImagePointer GetNormalVectorImage() const
    { return m_NormalVectorImage; }

  /** Get the lambda controlling the exponential decay used to calculate the
   *  weight value w from the distance to the closest border point.  Must be
   *  negative. */
  void SetLambda( WeightType l )
    { m_lambda = l; }
  const WeightType GetLambda() const
    { return m_lambda; }

  /** Get the weighting image */
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

  /** Get the normal and tangential components of the deformation field */
  virtual const OutputImagePointer GetOutputTangentialImage() const
    { return m_OutputTangentialImage; }
  virtual const OutputImagePointer GetOutputNormalImage() const
    { return m_OutputNormalImage; }

  /** Whether to compute the motion field regularization term (for testing)
   *  Default: true
   */
  void SetComputeRegularizationTerm( bool compute );
  bool GetComputeRegularizationTerm() const;

  /** Whether to compute the intensity distance term (for testing)
   *  Default: true
   */
  void SetComputeIntensityDistanceTerm( bool compute );
  bool GetComputeIntensityDistanceTerm() const;

  /** Whether to use diffusive regularization (if false, w=0 and uses typical
   *  Gaussian smoothing.  Default: true
   */
  void SetUseDiffusiveRegularization( bool diffuse )
    { m_UseDiffusiveRegularization = diffuse; }
  bool GetUseDiffusiveRegularization() const
    { return m_UseDiffusiveRegularization; }

protected:
  ImageToImageDiffusiveDeformableRegistrationFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Initialize the state of the filter and equation before each iteration. */
  virtual void InitializeIteration();

  /** This method populates an update buffer with changes for each pixel in the
   * output using the ThreadedCalculateChange() method and a multithreading
   * mechanism. Returns value is a time step to be used for the update. */
  virtual TimeStepType CalculateChange();

  /** Allocate the update buffer - reimplented here to also allocate other
   *  internal images. */
  virtual void AllocateUpdateBuffer();

  /** All other initialization before the registration loop */
  virtual void Initialize();

  /** Helper function to allocate an image based on a template */
  template< class UnallocatedImageType, class TemplateImageType >
  void AllocateSpaceForImage( UnallocatedImageType& inputImage,
                              const TemplateImageType& templateImage );

  /** Compute the normal vector image and weighting factor w given the
   *  surface border polydata.
   */
  virtual void ComputeNormalVectorAndWeightImages();

  /** Computes the weighting factor w from the distance to the border.  The
   *  weight should be 1 near the border and 0 away from the border, and
   *  weights between the diffusive regularization and a more typical Gaussian
   *  regularization.
   */
  virtual WeightType ComputeWeightFromDistance( WeightType distance );

  /** Compute the diffusion tensor image */
  virtual void ComputeDiffusionTensorImage();

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
  typename UpdateBufferType::Pointer    m_UpdateBuffer;

  /** Whether or not to use the diffusive regularization. */
  bool                                  m_UseDiffusiveRegularization;

  /** The border surface and derived image of normal vectors. */
  BorderSurfacePointer                  m_BorderSurface;
  BorderSurfacePointer                  m_BorderNormalsSurface;
  NormalVectorImagePointer              m_NormalVectorImage;

  /** The weighting image between the diffusive and gaussian regularizations */
  WeightImagePointer                            m_WeightImage;

  /** The lambda factor for computing the weight from distance.  Weight is
    * modeled as exponential decay: weight = e^(lambda * distance).
    * (lamba must be negative!)
    */
  WeightType                                    m_lambda;

  /** The image of the tangential diffusion tensors
   * Calculate div( P^TP gradient(motion field tangential component image) )
   */
  DiffusionTensorImagePointer           m_TangentialDiffusionTensorImage;

  /** The image of the normal diffusion tensors
   * Calculate div( w(n^T gradient(motion field normal component image) ) n )
   */
  DiffusionTensorImagePointer           m_NormalDiffusionTensorImage;

  /** Extracts the tangential and normal components of the deformation field */
  OutputImagePointer                    m_OutputTangentialImage;
  OutputImagePointer                    m_OutputNormalImage;

  /** Extracts the border surface normals from the border surface */
  typedef vtkPolyDataNormals                    PolyDataNormalsType;
  typedef vtkSmartPointer< vtkPolyDataNormals > PolyDataNormalsPointer;
  PolyDataNormalsPointer                        m_PolyDataNormals;

  /** Computes the nearest point on the polydata to a given point */
  typedef vtkPointLocator                       PointLocatorType;
  typedef PointLocatorType *                    PointLocatorPointer;
  PointLocatorPointer                           m_PointLocator;

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

/** TODO LIST - essential

  // Calculate w in ComputeNormalVectorAndWeightImages() - but don't bother if
  // m_UseDiffusiveRegularization == false (do it in the loop with the normals)

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
