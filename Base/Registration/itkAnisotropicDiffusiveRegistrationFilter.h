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
#ifndef __itkAnisotropicDiffusiveRegistrationFilter_h
#define __itkAnisotropicDiffusiveRegistrationFilter_h

#include "itkDiffusiveRegistrationFilter.h"

#include "vtkPolyData.h"
#include "vtkSmartPointer.h"

namespace itk
{

/** \class itkAnisotropicDiffusiveRegistrationFilter
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
class ITK_EXPORT AnisotropicDiffusiveRegistrationFilter
  : public DiffusiveRegistrationFilter< TFixedImage,
                                        TMovingImage,
                                        TDeformationField >
{
public:
  /** Standard class typedefs. */
  typedef AnisotropicDiffusiveRegistrationFilter            Self;
  typedef DiffusiveRegistrationFilter< TFixedImage,
                                       TMovingImage,
                                       TDeformationField >  Superclass;
  typedef SmartPointer< Self >                              Pointer;
  typedef SmartPointer< const Self >                        ConstPointer;

  /** Method for creation through the object factory. */
  itkDiffusiveRegistrationFilterNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, DiffusiveRegistrationFilter);

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
  typedef typename Superclass::ThreadRegionType         ThreadRegionType;

  /** Output image and update buffer types */
  typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;
  typedef typename Superclass::NeighborhoodType         NeighborhoodType;
  typedef typename Superclass::UpdateBufferRegionType   UpdateBufferRegionType;

  /** The registration function type */
  typedef AnisotropicDiffusiveRegistrationFunction
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

  /** Normal vector types */
  typedef typename RegistrationFunctionType::NormalVectorType
      NormalVectorType;
  typedef typename RegistrationFunctionType::NormalVectorImageType
      NormalVectorImageType;
  typedef typename NormalVectorImageType::Pointer
      NormalVectorImagePointer;
  typedef typename RegistrationFunctionType::NormalVectorNeighborhoodType
      NormalVectorNeighborhoodType;
  typedef itk::ImageRegionIterator< NormalVectorImageType >
      NormalVectorImageRegionType;
  typedef typename NormalVectorImageType::RegionType
      ThreadNormalVectorImageRegionType;

  /** Types for weighting between the anisotropic and diffusive (Gaussian)
    * regularization */
  typedef double                                        WeightType;
  typedef itk::Image< WeightType, ImageDimension >      WeightImageType;
  typedef typename WeightImageType::Pointer             WeightImagePointer;
  typedef itk::ImageRegionIterator< WeightImageType >   WeightImageRegionType;

  /** Organ boundary surface types */
  typedef vtkSmartPointer< vtkPolyData >                BorderSurfacePointer;

  /** Set/get whether to use the anisotropic diffusive regularization.  If
   *  false, the weighting term w=0 and Gaussian regularization is used.
   *  Default: true */
  void SetUseAnisotropicRegularization( bool diffuse )
    { this->GetRegistrationFunctionPointer()->
      SetUseAnisotropicRegularization( diffuse ); }
  bool GetUseAnisotropicRegularization() const
    { return this->GetRegistrationFunctionPointer()->
      GetUseAnisotropicRegularization(); }

  /** Set/get the organ boundary polydata, which must be in the same space as
   *  the fixed image.  Border normals are computed on this polydata, so it
   *  may be changed over the course of the registration. */
  virtual void SetBorderSurface( BorderSurfacePointer border )
    { m_BorderSurface = border; }
  virtual const BorderSurfacePointer GetBorderSurface() const
    { return m_BorderSurface; }

  /** Set/get the lambda that controls the exponential decay used to calculate
   *  the weight value w as a function of the distance to the closest border
   *  point.  Must be negative. */
  void SetLambda( WeightType l )
    { if ( l < 0 ) { m_lambda = l; } }
  WeightType GetLambda() const
    { return m_lambda; }

  /** Set/get the image of the normal vectors.  Setting the normal vector
   * image overrides the border surface polydata if a border surface was
   * also supplied. */
  virtual void SetNormalVectorImage( NormalVectorImagePointer normalImage )
    { m_NormalVectorImage = normalImage; }
  virtual const NormalVectorImagePointer GetNormalVectorImage() const
    { return m_NormalVectorImage; }

  /** Set/get the weighting image.  Setting the weighting image overrides
   * the border surface polydata and lambda if a border surface was also
   * supplied.  */
  virtual void SetWeightImage( WeightImagePointer weightImage )
    { m_WeightImage = weightImage; }
  virtual const WeightImagePointer GetWeightImage() const
    { return m_WeightImage; }

  /** Get the image of the normal diffusion tensors */
  virtual const DiffusionTensorImagePointer GetNormalDiffusionTensorImage()
    const
    { return m_NormalDiffusionTensorImage; }

  /** Get the image of the normal diffusion tensor derivatives */
  virtual const TensorDerivativeImagePointer
      GetNormalDiffusionTensorDerivativeImage() const
    { return m_NormalDiffusionTensorDerivativeImage; }

  /** Get the normal components of the deformation field */
  virtual const OutputImagePointer GetNormalDeformationFieldImage() const
    { return m_NormalDeformationField; }

protected:
  AnisotropicDiffusiveRegistrationFilter();
  virtual ~AnisotropicDiffusiveRegistrationFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Allocate images used during the registration. */
  virtual void AllocateImages();

  /** Compute the normals for the border surface. */
  void ComputeBorderSurfaceNormals();

  /** Computes the normal vector image and weighting factors w given the
   *  surface border polydata. */
  virtual void ComputeNormalVectorAndWeightImages( bool computeNormals,
                                                   bool computeWeights );

  /** Computes the weighting factor w from the distance to the border.  The
   *  weight should be 1 near the border and 0 away from the border. */
  virtual WeightType ComputeWeightFromDistance( WeightType distance );

  /** Computes the diffusion tensor images */
  virtual void ComputeDiffusionTensorImages();

  /** Computes the first derivatives of the diffusion tensor images */
  virtual void ComputeDiffusionTensorDerivativeImages();

  /** Updates the deformation vector component images */
  virtual void UpdateDeformationVectorComponentImages();

  /** This method populates an update buffer with changes for each pixel in the
   * output, using the ThreadedCalculateChange() method and a multithreading
   * mechanism. Return value is a time step to be used for the update.
   * \sa ThreadedCalculateChange */
  virtual TimeStepType CalculateChange();

  /** Inherited from superclass - do not call this function!  Call the other
   *  ThreadedCalculateChange function instead */
  virtual TimeStepType ThreadedCalculateChange(
      const ThreadRegionType &regionToProcess,
      const ThreadDiffusionTensorImageRegionType &tensorRegionToProcess,
      const ThreadTensorDerivativeImageRegionType
        &tensorDerivativeRegionToProcess,
      const ThreadDeformationVectorComponentImageRegionType
        &deformationComponentRegionToProcess,
      int threadId );

  /** Does the actual work of calculating change over a region supplied by
   * the multithreading mechanism.
   * \sa CalculateChange
   * \sa CalculateChangeThreaderCallback */
  virtual TimeStepType ThreadedCalculateChange(
      const ThreadRegionType &regionToProcess,
      const ThreadNormalVectorImageRegionType &normalVectorRegionToProcess,
      const ThreadDiffusionTensorImageRegionType &tensorRegionToProcess,
      const ThreadTensorDerivativeImageRegionType
        &tensorDerivativeRegionToProcess,
      const ThreadDeformationVectorComponentImageRegionType
        &deformationComponentRegionToProcess,
      int threadId );

  /** Create the registration function, with default parameters for
    * ComputeRegularizationTerm and ComputeIntensityDistanceTerm. */
  virtual void CreateRegistrationFunction();

  /** Get the registration function pointer */
  virtual RegistrationFunctionType * GetRegistrationFunctionPointer() const;

private:
  // Purposely not implemented
  AnisotropicDiffusiveRegistrationFilter(const Self&);
  void operator=(const Self&); // Purposely not implemented

  /** Structure for passing information into static callback methods.  Used in
   * the subclasses threading mechanisms. */
  struct DenseFDThreadStruct
    {
    AnisotropicDiffusiveRegistrationFilter *Filter;
    TimeStepType TimeStep;
    TimeStepType *TimeStepList;
    bool *ValidTimeStepList;
    };

  /** This callback method uses SplitUpdateContainer to acquire a region
  * which it then passes to ThreadedCalculateChange for processing. */
  static ITK_THREAD_RETURN_TYPE CalculateChangeThreaderCallback( void *arg );

  /** Organ boundary surface and surface of border normals */
  BorderSurfacePointer                m_BorderSurface;

  /** Image storing information we will need for each voxel on every
   *  registration iteration */
  NormalVectorImagePointer            m_NormalVectorImage;
  WeightImagePointer                  m_WeightImage;
  OutputImagePointer                  m_NormalDeformationField;
  DiffusionTensorImagePointer         m_NormalDiffusionTensorImage;
  TensorDerivativeImagePointer        m_NormalDiffusionTensorDerivativeImage;
  DeformationComponentImageArrayType  m_NormalDeformationComponentImages;

  /** The lambda factor for computing the weight from distance.  Weight is
   * modeled as exponential decay: weight = e^(lambda * distance).
   * (lamba must be negative) */
  WeightType                          m_lambda;

};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkAnisotropicDiffusiveRegistrationFilter+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkAnisotropicDiffusiveRegistrationFilter.txx"
#endif

#endif
