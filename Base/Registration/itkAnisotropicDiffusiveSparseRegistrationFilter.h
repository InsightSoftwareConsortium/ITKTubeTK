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
#ifndef __itkAnisotropicDiffusiveSparseRegistrationFilter_h
#define __itkAnisotropicDiffusiveSparseRegistrationFilter_h

#include "itkDiffusiveRegistrationFilter.h"

namespace itk
{

/** \class itkAnisotropicDiffusiveSparseRegistrationFilter
 * \brief Registration filter for registrations using anisotropic diffusive
 * regularizers, for example for sliding organ registration.
 *
 * Traditional deformable image registration imposes a uniform
 * smoothness constraint on the deformation field. This
 * is not appropriate when registering images visualizing organs
 * that slide relative to each other, and therefore leads to registration
 * inaccuracies.
 *
 * This filter includes a regularization term based on anisotropic diffusion
 * that accomodates deformation field discontinuities that are expected when
 * considering sliding motion.
 *
 * TODO
 *
 * See: D.F. Pace et al., Deformable image registration of sliding organs using
 * anisotropic diffusive regularization, ISBI 2011.
 *
 * This class is templated over the type of the fixed image, the type of the
 * moving image and the type of the deformation field.
 *
 * \sa itkDiffusiveRegistrationFilter
 * \sa itkAnisotropicDiffusiveRegistrationFilter
 * \sa itkAnisotropicDiffusiveRegistrationFunction
 * \ingroup DeformableImageRegistration
 * \ingroup MultiThreaded
 */

template < class TFixedImage, class TMovingImage, class TDeformationField >
class ITK_EXPORT AnisotropicDiffusiveSparseRegistrationFilter
  : public DiffusiveRegistrationFilter< TFixedImage,
                                        TMovingImage,
                                        TDeformationField >
{
public:
  /** Standard class typedefs. */
  typedef AnisotropicDiffusiveSparseRegistrationFilter      Self;
  typedef DiffusiveRegistrationFilter< TFixedImage,
                                       TMovingImage,
                                       TDeformationField >  Superclass;
  typedef SmartPointer< Self >                              Pointer;
  typedef SmartPointer< const Self >                        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

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
  typedef typename Superclass::OutputImageType          OutputImageType;
  typedef typename Superclass::OutputImagePointer       OutputImagePointer;
  typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;

  /** The registration function type */
  typedef typename Superclass::RegistrationFunctionType
      RegistrationFunctionType;
  typedef typename Superclass::SpacingType              SpacingType;

  /** Deformation field types. */
  typedef typename Superclass::DeformationVectorType    DeformationVectorType;
  typedef typename Superclass::DeformationVectorComponentType
      DeformationVectorComponentType;
  typedef typename Superclass::DeformationComponentImageArrayType
      DeformationComponentImageArrayType;

  /** Diffusion tensor image types */
  typedef typename Superclass::DiffusionTensorType      DiffusionTensorType;
  typedef typename Superclass::DiffusionTensorImageType
      DiffusionTensorImageType;

  /** Scalar derivative image types */
  typedef typename Superclass::ScalarDerivativeImageType
      ScalarDerivativeImageType;
  typedef typename Superclass::ScalarDerivativeImagePointer
      ScalarDerivativeImagePointer;
  typedef typename Superclass::ScalarDerivativeImageArrayType
      ScalarDerivativeImageArrayType;

  /** Tensor derivative matrix image types */
  typedef typename Superclass::TensorDerivativeImageType
      TensorDerivativeImageType;
  typedef typename Superclass::TensorDerivativeImagePointer
      TensorDerivativeImagePointer;
  typedef typename Superclass::TensorDerivativeImageArrayType
      TensorDerivativeImageArrayType;

  /** Typedefs for the multiplication vectors */
  typedef typename Superclass::DeformationVectorImageArrayType
      DeformationVectorImageArrayType;
  typedef typename Superclass::DeformationVectorImageRegionType
      DeformationVectorImageRegionType;
  typedef typename
      RegistrationFunctionType::DeformationVectorImageRegionArrayType
      DeformationVectorImageRegionArrayType;

  /** Normal vector types.  There are three normals at each voxel, which are
   *  stored in a matrix.  If the normals are based on the structure tensor,
   *  then the matrix will be symmetric, but we won't enforce that. */
  typedef double NormalVectorComponentType;
  typedef typename itk::Vector< NormalVectorComponentType, ImageDimension >
      NormalVectorType;
  typedef typename itk::Matrix< NormalVectorComponentType,
                                ImageDimension,
                                ImageDimension >          NormalMatrixType;
  typedef itk::Image< NormalMatrixType, ImageDimension >
      NormalMatrixImageType;
  typedef typename NormalMatrixImageType::Pointer
      NormalMatrixImagePointer;
  typedef ZeroFluxNeumannBoundaryCondition< NormalMatrixImageType >
      NormalMatrixImageBoundaryConditionType;
  typedef itk::ImageRegionIterator< NormalMatrixImageType >
      NormalMatrixImageRegionType;

  /** Types for weighting between the anisotropic and diffusive (Gaussian)
    * regularization - also a matrix, likely symmetric but we won't enforce
    * that here. */
  typedef double                                        WeightComponentType;
  typedef typename itk::Matrix< WeightComponentType,
                                ImageDimension,
                                ImageDimension >        WeightType;
  typedef itk::Image< WeightType, ImageDimension >      WeightImageType;
  typedef typename WeightImageType::Pointer             WeightImagePointer;
  typedef itk::ImageRegionIterator< WeightImageType >   WeightImageRegionType;

  /** Organ boundary surface types */
  typedef vtkPolyData                                   BorderSurfaceType;
  typedef vtkSmartPointer< BorderSurfaceType >          BorderSurfacePointer;

  /** The number of div(Tensor \grad u)v terms we sum for the regularizer.
   *  Reimplement in derived classes. */
  virtual int GetNumberOfTerms() const
    { return 4; }

  /** Set/get the organ boundary polydata, which must be in the same space as
   *  the fixed image.  Border normals are computed on this polydata, so it
   *  may be changed over the course of the registration. */
  virtual void SetBorderSurface( BorderSurfaceType * border )
    { m_BorderSurface = border; }
  virtual BorderSurfaceType * GetBorderSurface() const
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
  virtual void SetNormalMatrixImage( NormalMatrixImageType * normalImage )
    { m_NormalMatrixImage = normalImage; }
  virtual const NormalMatrixImageType * GetNormalMatrixImage() const
    { return m_NormalMatrixImage; }

  /** Set/get the weighting image.  Setting the weighting image overrides
   * the border surface polydata and lambda if a border surface was also
   * supplied.  */
  virtual void SetWeightImage( WeightImageType * weightImage )
    { m_WeightImage = weightImage; }
  virtual const WeightImageType * GetWeightImage() const
    { return m_WeightImage; }

  /** Get the normal components of the deformation field.  The normal
   *  deformation field component images are the same for both the SMOOTH_NORMAL
   *  and PROP_NORMAL terms, so we will return one arbitrarily. */
  virtual const DeformationFieldType * GetNormalDeformationComponentImage()
      const
    {
    return this->GetDeformationComponentImage( SMOOTH_NORMAL );
    }

protected:
  AnisotropicDiffusiveSparseRegistrationFilter();
  virtual ~AnisotropicDiffusiveSparseRegistrationFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Handy for array indexing. */
  enum DivTerm { SMOOTH_TANGENTIAL,
                 SMOOTH_NORMAL,
                 PROP_TANGENTIAL,
                 PROP_NORMAL };

  /** Allocate the deformation component images and their derivative images.
   *  (which may be updated throughout the registration). Reimplement in derived
   *  classes. */
  virtual void InitializeDeformationComponentAndDerivativeImages();

  /** Allocate and populate the diffusion tensor images.
   *  Reimplement in derived classes. */
  virtual void ComputeDiffusionTensorImages();

  /** Allocate and populate the images of multiplication vectors that the
   *  div(T \grad(u)) values are multiplied by.  Allocate and populate all or
   *  some of the multiplication vector images in derived classes.  Otherwise,
   *  default to e_l, where e_l is the lth canonical unit vector. */
  virtual void ComputeMultiplicationVectorImages();

  /** Updates the deformation vector component images on each iteration. */
  virtual void UpdateDeformationComponentImages();

  /** Computes the first- and second-order partial derivatives of the
   *  deformation component images on each iteration.  Override in derived
   *  classes if the deformation components image pointers are not unique, to
   *  avoid computing the same derivatives multiple times. */
  virtual void ComputeDeformationComponentDerivativeImages();

  /** If needed, allocates and computes the normal vector and weight images. */
  virtual void SetupNormalMatrixAndWeightImages();

  /** Compute the normals for the border surface. */
  void ComputeBorderSurfaceNormals();

  /** Computes the normal vector image and weighting factors w given the
   *  surface border polydata. */
  virtual void ComputeNormalMatrixAndWeightImages( bool computeNormals,
                                                   bool computeWeights );

  /** Computes the weighting factor w from the distance to the border.  The
   *  weight should be 1 near the border and 0 away from the border. */
  virtual WeightComponentType ComputeWeightFromDistance(
      const WeightComponentType distance ) const;

private:
  // Purposely not implemented
  AnisotropicDiffusiveSparseRegistrationFilter(const Self&);
  void operator=(const Self&); // Purposely not implemented

  /** Organ boundary surface and surface of border normals */
  BorderSurfacePointer                m_BorderSurface;

  /** Image storing information we will need for each voxel on every
   *  registration iteration */
  NormalMatrixImagePointer            m_NormalMatrixImage;
  WeightImagePointer                  m_WeightImage;

  /** The lambda factor for computing the weight from distance.  Weight is
   * modeled as exponential decay: weight = e^(lambda * distance).
   * (lamba must be negative) */
  WeightComponentType                 m_lambda;

};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkAnisotropicDiffusiveSparseRegistrationFilter+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkAnisotropicDiffusiveSparseRegistrationFilter.txx"
#endif

#endif
