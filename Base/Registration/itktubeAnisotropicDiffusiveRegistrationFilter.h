/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeAnisotropicDiffusiveRegistrationFilter_h
#define __itktubeAnisotropicDiffusiveRegistrationFilter_h

#include "itktubeDiffusiveRegistrationFilter.h"

#include <vtkSmartPointer.h>

class vtkFloatArray;
class vtkPointLocator;
class vtkPolyData;

namespace itk
{

namespace tube
{

/** \class AnisotropicDiffusiveRegistrationFilter
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
 * that accommodates deformation field discontinuities that are expected when
 * considering sliding motion.
 *
 * The regularization term uses a specified border between the organs
 * ( stored as a vtkPolyData * ) and enforces coupling between the organs while
 * allowing the motion field to exhibit sliding motion at the organ interface.
 *
 * See: D.F. Pace et al., Deformable image registration of sliding organs using
 * anisotropic diffusive regularization, ISBI 2011.
 *
 * This class is templated over the type of the fixed image, the type of the
 * moving image and the type of the deformation field.
 *
 * \sa DiffusiveRegistrationFilter
 * \sa AnisotropicDiffusiveRegistrationFunction
 * \ingroup DeformableImageRegistration
 * \ingroup MultiThreaded
 */

template< class TFixedImage, class TMovingImage, class TDeformationField >
class AnisotropicDiffusiveRegistrationFilter
  : public DiffusiveRegistrationFilter< TFixedImage, TMovingImage,
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
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( AnisotropicDiffusiveRegistrationFilter,
                DiffusiveRegistrationFilter );

  /** Inherit some parameters from the superclass. */
  itkStaticConstMacro( ImageDimension, unsigned int, Superclass::ImageDimension );

  /** Convenient typedefs from the superclass. */
  typedef typename Superclass::FixedImageType
      FixedImageType;
  typedef typename Superclass::FixedImagePointer
      FixedImagePointer;
  typedef typename Superclass::MovingImageType
      MovingImageType;
  typedef typename Superclass::MovingImagePointer
      MovingImagePointer;
  typedef typename Superclass::DeformationFieldType
      DeformationFieldType;
  typedef typename Superclass::DeformationFieldPointer
      DeformationFieldPointer;
  typedef typename Superclass::TimeStepType
      TimeStepType;
  typedef typename Superclass::FiniteDifferenceFunctionType
      FiniteDifferenceFunctionType;
  typedef typename Superclass::OutputImageType
      OutputImageType;
  typedef typename Superclass::OutputImagePointer
      OutputImagePointer;
  typedef typename Superclass::OutputImageRegionType
      OutputImageRegionType;

  /** The registration function type */
  typedef typename Superclass::RegistrationFunctionType
      RegistrationFunctionType;

  /** Deformation field types. */
  typedef typename Superclass::DeformationVectorType
      DeformationVectorType;
  typedef typename Superclass::DeformationVectorComponentType
      DeformationVectorComponentType;

  /** Diffusion tensor image types */
  typedef typename Superclass::DiffusionTensorType
      DiffusionTensorType;
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
  typedef itk::ImageRegionIterator< NormalVectorImageType >
      NormalVectorImageRegionType;
  typedef typename NormalVectorImageType::RegionType
      ThreadNormalVectorImageRegionType;

  /** Types for weighting between the anisotropic and diffusive ( Gaussian )
    * regularization */
  typedef double                                        WeightType;
  typedef itk::Image< WeightType, ImageDimension >      WeightImageType;
  typedef typename WeightImageType::Pointer             WeightImagePointer;
  typedef itk::ImageRegionIterator< WeightImageType >   WeightImageRegionType;
  typedef typename WeightImageType::RegionType
      ThreadWeightImageRegionType;

  /** Organ boundary surface types */
  typedef vtkPolyData                                   BorderSurfaceType;
  typedef vtkSmartPointer< BorderSurfaceType >          BorderSurfacePointer;

  /** The number of div( Tensor \grad u )v terms we sum for the regularizer.
   *  Reimplement in derived classes. */
  virtual int GetNumberOfTerms( void ) const
    { return 2; }

  /** Set/get the organ boundary polydata, which must be in the same space as
   *  the fixed image.  Border normals are computed on this polydata, so it
   *  may be changed over the course of the registration. */
  virtual void SetBorderSurface( BorderSurfaceType * border )
    { m_BorderSurface = border; }
  virtual BorderSurfaceType * GetBorderSurface( void ) const
    { return m_BorderSurface; }

  /** Set/get the lambda that controls the decay of the weight value w as a
   *  function of the distance to the closest border point.  If gamma=-1, then
   *  w decays exponentially ( w = e^( -1.0*lambda*distance ) ).  Otherwise, w
   *  decays exponentially using a Dirac-shaped function
   *  ( w = 1 / ( 1 + lambda*gamma*e^( -1.0*lambda*distance^2 ) ) ).  Lambda must
   *  be positive. */
  void SetLambda( WeightType l )
    { if( l > 0 ) { m_Lambda = l; } }
  WeightType GetLambda( void ) const
    { return m_Lambda; }

  /** Set/get the gamma that controls the decay of the weight value w as a
   *  function of the distance to the closest border point.  If gamma=-1, then
   *  w decays exponentially ( w = e^( -1.0*lambda*distance ) ).  Otherwise, w
   *  decays exponentially using a Dirac-shaped function
   *  ( w = 1 / ( 1 + lambda*gamma*e^( -1.0*lambda*distance^2 ) ) ).  Gamma must
   *  be positive or -1.0. */
  void SetGamma( WeightType g )
    { if( g > 0 || g == -1.0 ) { m_Gamma = g; } }
  WeightType GetGamma( void ) const
    { return m_Gamma; }

  /** Set/get the image of the normal vectors.  Setting the normal vector
   * image overrides the border surface polydata if a border surface was
   * also supplied. */
  virtual void SetNormalVectorImage( NormalVectorImageType * normalImage )
    { m_NormalVectorImage = normalImage; }
  virtual NormalVectorImageType * GetNormalVectorImage( void ) const
    { return m_NormalVectorImage; }
  virtual NormalVectorImageType * GetHighResolutionNormalVectorImage( void ) const
    { return m_HighResolutionNormalVectorImage; }

  /** Set/get the weighting image.  Setting the weighting image overrides
   * the border surface polydata and lambda/gamma if a border surface was also
   * supplied.  */
  virtual void SetWeightImage( WeightImageType * weightImage )
    { m_WeightImage = weightImage; }
  virtual WeightImageType * GetWeightImage( void ) const
    { return m_WeightImage; }
  virtual WeightImageType * GetHighResolutionWeightImage( void ) const
    { return m_HighResolutionWeightImage; }

  /** Get the normal components of the deformation field. */
  virtual const DeformationFieldType * GetNormalDeformationComponentImage( void )
      const
    {
    return this->GetDeformationComponentImage( NORMAL );
    }

protected:
  AnisotropicDiffusiveRegistrationFilter( void );
  virtual ~AnisotropicDiffusiveRegistrationFilter( void ) {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

  /** Handy for array indexing. */
  enum DivTerm { TANGENTIAL, NORMAL };

  /** Allocate the deformation component images and their derivative images.
   *  ( which may be updated throughout the registration ). Reimplement in derived
   *  classes. */
  virtual void InitializeDeformationComponentAndDerivativeImages( void );

  /** Allocate and populate the diffusion tensor images.
   *  Reimplement in derived classes. */
  virtual void ComputeDiffusionTensorImages( void );

  /** Allocate and populate the images of multiplication vectors that the
   *  div( T \grad( u ) ) values are multiplied by.  Allocate and populate all or
   *  some of the multiplication vector images in derived classes.  Otherwise,
   *  default to e_l, where e_l is the lth canonical unit vector. */
  virtual void ComputeMultiplicationVectorImages( void );

  /** Updates the deformation vector component images on each iteration. */
  virtual void UpdateDeformationComponentImages( OutputImageType * output );

  /** If needed, allocates and computes the normal vector and weight images. */
  virtual void SetupNormalVectorAndWeightImages( void );

  /** Compute the normals for the border surface. */
  void ComputeBorderSurfaceNormals( void );

  /** Computes the normal vector image and weighting factors w given the
   *  surface border polydata. */
  virtual void ComputeNormalVectorAndWeightImages( bool computeNormals,
                                                   bool computeWeights );

  /** Computes the normal vectors and distances to the closest point given
   *  an initialized vtkPointLocator and the surface border normals */
  virtual void GetNormalsAndDistancesFromClosestSurfacePoint( 
      bool computeNormals, bool computeWeights );

  /** Does the actual work of updating the output over an output region supplied
   *  by the multithreading mechanism.
   *  \sa GetNormalsAndDistancesFromClosestSurfacePoint
   *  \sa GetNormalsAndDistancesFromClosestSurfacePointThreaderCallback */
  virtual void ThreadedGetNormalsAndDistancesFromClosestSurfacePoint( 
      vtkPointLocator * pointLocator,
      vtkFloatArray * normalData,
      ThreadNormalVectorImageRegionType & normalRegionToProcess,
      ThreadWeightImageRegionType & weightRegionToProcess,
      bool computeNormals,
      bool computeWeights,
      int threadId );

  /** Computes the weighting factor w from the distance to the border using
   *  exponential decay.  The weight should be 1 near the border and 0 away from
   *  the border. */
  virtual WeightType ComputeWeightFromDistanceExponential( 
      const WeightType distance ) const;

  /** Computes the weighting factor w from the distance to the border using
   *  a Dirac-shaped function.  The weight should be 1 near the border and 0
   *  away from the border. */
  virtual WeightType ComputeWeightFromDistanceDirac( 
      const WeightType distance ) const;

private:
  // Purposely not implemented
  AnisotropicDiffusiveRegistrationFilter( const Self& );
  void operator=( const Self& ); // Purposely not implemented

  /** Structure for passing information into static callback methods.  Used in
   * the subclasses threading mechanisms. */
  struct AnisotropicDiffusiveRegistrationFilterThreadStruct
    {
    AnisotropicDiffusiveRegistrationFilter * Filter;
    vtkPointLocator * PointLocator;
    vtkFloatArray * NormalData;
    ThreadNormalVectorImageRegionType NormalVectorImageLargestPossibleRegion;
    ThreadWeightImageRegionType WeightImageLargestPossibleRegion;
    bool ComputeNormals;
    bool ComputeWeights;
    }; // End struct AnisotropicDiffusiveRegistrationFilterThreadStruct

  /** This callback method uses ImageSource::SplitRequestedRegion to acquire an
   * output region that it passes to
   * ThreadedGetNormalsAndDistancesFromClosestSurfacePoint for processing. */
  static ITK_THREAD_RETURN_TYPE
      GetNormalsAndDistancesFromClosestSurfacePointThreaderCallback( 
          void * arg );

  /** Organ boundary surface and surface of border normals */
  BorderSurfacePointer                m_BorderSurface;

  /** Image storing information we will need for each voxel on every
   *  registration iteration */
  NormalVectorImagePointer            m_NormalVectorImage;
  WeightImagePointer                  m_WeightImage;

  /** Highest resolution versions of the normal and weight images, useful
   *  to calculate once ( setting m_ImageAttributeImage ) at the highest
   *  resolution during multiresolution registration, and then resampling on
   *  each scale.  The normal matrix image and weight structures image are
   *  resampled using nearest neighbor, while the weight regularizations image
   *  are resampled using a linear interpolation */
  NormalVectorImagePointer            m_HighResolutionNormalVectorImage;
  WeightImagePointer                  m_HighResolutionWeightImage;

  /** The lambda/gamma factors for computing the weight from distance. */
  WeightType                          m_Lambda;
  WeightType                          m_Gamma;

}; // End class AnisotropicDiffusiveRegistrationFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeAnisotropicDiffusiveRegistrationFilter.hxx"
#endif

#endif // End !defined( __itktubeAnisotropicDiffusiveRegistrationFilter_h )
