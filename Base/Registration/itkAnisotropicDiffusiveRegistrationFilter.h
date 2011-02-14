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

  /** Typedefs used in multithreading */
  typedef typename Superclass::OutputImageType          OutputImageType;
  typedef typename Superclass::OutputImagePointer       OutputImagePointer;
  typedef typename Superclass::UpdateBufferType         UpdateBufferType;

  /** Output image and update buffer types */
  typedef typename Superclass::OutputImageRegionType    OutputImageRegionType;
  typedef typename Superclass::NeighborhoodType         NeighborhoodType;
  typedef typename Superclass::UpdateBufferRegionType   UpdateBufferRegionType;

  /** The registration function type */
  typedef typename Superclass::RegistrationFunctionType
      RegistrationFunctionType;
  typedef typename Superclass::RegistrationFunctionPointer
      RegistrationFunctionPointer;
  typedef typename Superclass::RegularizationFunctionPointer
      RegularizationFunctionPointer;
  typedef typename Superclass::SpacingType              SpacingType;

  /** Deformation field types. */
  typedef typename Superclass::DeformationVectorType
      DeformationVectorType;
  typedef typename Superclass::DeformationVectorComponentType
      DeformationVectorComponentType;
  typedef typename Superclass::DeformationVectorComponentImageType
      DeformationVectorComponentImageType;
  typedef typename Superclass::DeformationVectorComponentImagePointer
      DeformationVectorComponentImagePointer;
  typedef typename Superclass::DeformationComponentImageArrayType
      DeformationComponentImageArrayType;
  typedef typename Superclass::DeformationVectorComponentNeighborhoodType
      DeformationVectorComponentNeighborhoodType;
  typedef typename Superclass::DeformationVectorComponentNeighborhoodArrayType
      DeformationVectorComponentNeighborhoodArrayType;

  /** Diffusion tensor image types */
  typedef typename Superclass::DiffusionTensorImageType
      DiffusionTensorImageType;
  typedef typename Superclass::DiffusionTensorImagePointer
      DiffusionTensorImagePointer;
  typedef typename Superclass::DiffusionTensorNeighborhoodType
      DiffusionTensorNeighborhoodType;

  /** Tensor derivative matrix image types */
  typedef typename Superclass::TensorDerivativeImageType
      TensorDerivativeImageType;
  typedef typename Superclass::TensorDerivativeImagePointer
      TensorDerivativeImagePointer;
  typedef typename Superclass::TensorDerivativeImageRegionType
      TensorDerivativeImageRegionType;

  /** Typedefs for the multiplication vectors */
  typedef typename RegistrationFunctionType::DeformationVectorImageRegionType
      DeformationVectorImageRegionType;
  typedef typename
      RegistrationFunctionType::DeformationVectorImageRegionArrayType
      DeformationVectorImageRegionArrayType;
  typedef typename
      RegistrationFunctionType::DeformationVectorImageRegionArrayArrayType
      DeformationVectorImageRegionArrayArrayType;

  /** Normal vector types */
  typedef double NormalVectorComponentType;
  typedef itk::Vector< NormalVectorComponentType, ImageDimension >
      NormalVectorType;
  typedef itk::Image< NormalVectorType, ImageDimension >
      NormalVectorImageType;
  typedef typename NormalVectorImageType::Pointer
      NormalVectorImagePointer;
  typedef ZeroFluxNeumannBoundaryCondition< NormalVectorImageType >
      NormalVectorImageBoundaryConditionType;
  typedef ConstNeighborhoodIterator
      < NormalVectorImageType, NormalVectorImageBoundaryConditionType >
      NormalVectorNeighborhoodType;
  typedef itk::ImageRegionIterator< NormalVectorImageType >
      NormalVectorImageRegionType;

  /** Types for weighting between the anisotropic and diffusive (Gaussian)
    * regularization */
  typedef double                                        WeightType;
  typedef itk::Image< WeightType, ImageDimension >      WeightImageType;
  typedef typename WeightImageType::Pointer             WeightImagePointer;
  typedef itk::ImageRegionIterator< WeightImageType >   WeightImageRegionType;

  /** Organ boundary surface types */
  typedef vtkPolyData                                   BorderSurfaceType;
  typedef vtkSmartPointer< BorderSurfaceType >          BorderSurfacePointer;

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
  virtual void SetNormalVectorImage( NormalVectorImageType * normalImage )
    { m_NormalVectorImage = normalImage; }
  virtual const NormalVectorImageType * GetNormalVectorImage() const
    { return m_NormalVectorImage; }

  /** Set/get the weighting image.  Setting the weighting image overrides
   * the border surface polydata and lambda if a border surface was also
   * supplied.  */
  virtual void SetWeightImage( WeightImageType * weightImage )
    { m_WeightImage = weightImage; }
  virtual const WeightImageType * GetWeightImage() const
    { return m_WeightImage; }

  /** The number of div(Tensor \grad u)v terms we sum for the regularizer.
   *  Reimplement in derived classes. */
  virtual int GetNumberOfTerms() const
    { return 2; }

  /** Get the normal components of the deformation field. */
  virtual DeformationFieldType * GetNormalDeformationFieldImage() const
    {
    return this->GetDeformationFieldComponentImage( NORMAL );
    }

protected:
  AnisotropicDiffusiveRegistrationFilter();
  virtual ~AnisotropicDiffusiveRegistrationFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  enum DivTerm { TANGENTIAL, NORMAL };

  /** Initialize images used during the registration. */
  virtual void InitializeImageArrays();

  /** Compute the normals for the border surface. */
  void ComputeBorderSurfaceNormals();

  /** Computes the normal vector image and weighting factors w given the
   *  surface border polydata. */
  virtual void InitializeNormalVectorAndWeightImages( bool computeNormals,
                                                      bool computeWeights );

  /** Computes the weighting factor w from the distance to the border.  The
   *  weight should be 1 near the border and 0 away from the border. */
  virtual WeightType ComputeWeightFromDistance( const WeightType distance )
      const;

  /** Computes the diffusion tensor images */
  virtual void ComputeDiffusionTensorImages();

  /** Updates the deformation field component images */
  virtual void UpdateDeformationFieldComponentImages();

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

  /** Organ boundary surface and surface of border normals */
  BorderSurfacePointer                m_BorderSurface;

  /** Image storing information we will need for each voxel on every
   *  registration iteration */
  NormalVectorImagePointer            m_NormalVectorImage;
  WeightImagePointer                  m_WeightImage;

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
