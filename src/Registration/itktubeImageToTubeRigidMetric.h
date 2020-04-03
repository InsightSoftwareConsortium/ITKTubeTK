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

#ifndef __itktubeImageToTubeRigidMetric_h
#define __itktubeImageToTubeRigidMetric_h

#include <itkEuler3DTransform.h>
#include <itkCompensatedSummation.h>
#include <itkGaussianDerivativeImageFunction.h>
#include <itkImageToSpatialObjectMetric.h>

namespace itk
{

namespace tube
{

/** \class ImageToTubeRigidMetric
 * \brief Computes similarity between two objects to be registered
 * The metric implemented here corresponds to the following paper:
 * \link
 * http://www.cs.unc.edu/Research/MIDAG/pubs/papers/MICCAI01-aylwardVReg.pdf
 * The metric is based on the fact that vessel centerlines are scaled
 * intensity ridges in the image.
 *
 * \tparam TFixedImage Type of the Image to register against.
 * \tparam TMovingSpatialObject Type of the SpatialObject to register with,
 * could be a Tube, Group, etc.
 * \tparam TTubeSpatialObject Type of the tubes contained within the input
 * TMovingSpatialObject to use for the registration.
 *
 * \warning ( Derivative )
 */

template< class TFixedImage,
          class TMovingSpatialObject,
          class TTubeSpatialObject >
class ImageToTubeRigidMetric
  : public ImageToSpatialObjectMetric< TFixedImage, TMovingSpatialObject >
{
public:
  /** Standard class typedefs. */
  typedef ImageToTubeRigidMetric                Self;
  typedef ImageToSpatialObjectMetric< TFixedImage, TMovingSpatialObject >
                                                Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;

  /**  Dimension of the image and tube.  */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TFixedImage::ImageDimension );
  itkStaticConstMacro( TubeDimension, unsigned int,
    TTubeSpatialObject::ObjectDimension );

  typedef TFixedImage                           FixedImageType;
  typedef TMovingSpatialObject                  TubeTreeType;
  typedef TTubeSpatialObject                    TubeType;
  typedef typename TubeType::TubePointType      TubePointType;

  typedef double                                ScalarType;
  typedef GaussianDerivativeImageFunction< TFixedImage >
                                                DerivativeImageFunctionType;
  typedef typename Superclass::DerivativeType   DerivativeType;
  typedef typename Superclass::ParametersType   ParametersType;
  typedef typename Superclass::MeasureType      MeasureType;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ImageToTubeRigidMetric, ImageToSpatialObjectMetric );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  inline unsigned int GetNumberOfParameters( void ) const override
    {
    return this->m_Transform->GetNumberOfParameters();
    }

  /** Type used for representing point components  */
  typedef typename Superclass::CoordinateRepresentationType
                                               CoordinateRepresentationType;

  /** Type definition for the size */
  typedef typename TFixedImage::SizeType       SizeType;

  /** Type definition for the pixel type */
  typedef typename TFixedImage::PixelType      PixelType;

  /**  Type of the Transform Base class */
  typedef Euler3DTransform< ScalarType >           TransformType;
  typedef typename TransformType::Pointer          TransformPointer;
  typedef typename TransformType::InputPointType   InputPointType;
  typedef typename TransformType::OutputPointType  OutputPointType;
  typedef typename TransformType::ParametersType   TransformParametersType;
  typedef typename TransformType::JacobianType     TransformJacobianType;
  typedef TransformType::ParametersType            FeatureWeightsType;

  /** Get the Derivatives of the Match Measure */
  const DerivativeType & GetDerivative( const ParametersType &
    parameters ) const;
  void GetDerivative( const ParametersType & parameters,
    DerivativeType & derivative ) const override;

  /** Get the Value for SingleValue Optimizers */
  MeasureType  GetValue( const ParametersType & parameters ) const override;

  /** Get Value and Derivatives for MultipleValuedOptimizers */
  void GetValueAndDerivative( const ParametersType & parameters,
    MeasureType & Value, DerivativeType  & Derivative ) const override;

  /** Initialize the metric */
  void Initialize( void ) override;

  /** Control the radius scaling of the metric. */
  itkSetMacro( Kappa, ScalarType );
  itkGetConstMacro( Kappa, ScalarType );

  /** Set/Get the minimum scaling radius. */
  itkSetMacro( MinimumScalingRadius, ScalarType );
  itkGetConstMacro( MinimumScalingRadius, ScalarType );

  /** Set/Get the extent of the blurring calculation given in Gaussian
   * sigma's. */
  itkSetMacro( Extent, ScalarType );
  itkGetConstMacro( Extent, ScalarType );

  /** Set/Get the scalar weights associated with every point in the tube.
   * The index of the point weights should correspond to "standard tube tree
   * interation". */
  void SetFeatureWeights( FeatureWeightsType & featureWeights );
  itkGetConstReferenceMacro( FeatureWeights, FeatureWeightsType )

  TransformPointer GetTransform( void ) const
    {
    return dynamic_cast<TransformType*>( this->m_Transform.GetPointer() );
    }

  /** Downsample the tube points by this integer value. */

protected:
  ImageToTubeRigidMetric( void );
  virtual ~ImageToTubeRigidMetric( void );

  typedef CovariantVector< ScalarType, TubeDimension >    CovariantVectorType;
  typedef Vector< ScalarType, TubeDimension >             VectorType;

  typedef Matrix< ScalarType, TubeDimension, TubeDimension >   MatrixType; 

  typedef typename TubePointType::PointType        PointType;
  typedef vnl_vector< ScalarType >                 VnlVectorType;
  typedef vnl_matrix< ScalarType >                 VnlMatrixType;
  typedef CompensatedSummation< ScalarType >       CompensatedSummationType;

  virtual void ComputeCenterOfRotation( void );
  SizeValueType CountTubePoints( void );

  void GetDeltaAngles( const OutputPointType & x,
    const VnlVectorType & dx,
    const VectorType & offsets,
    ScalarType angle[3] ) const;

private:
  ImageToTubeRigidMetric( const Self& ); // purposely not implemented
  void operator=( const Self& ); // purposely not implemented

  typename DerivativeImageFunctionType::Pointer m_DerivativeImageFunction;

  ScalarType m_Kappa;
  ScalarType m_MinimumScalingRadius;
  ScalarType m_Extent;

  /** The center of rotation of the weighted tube points. */
  typedef PointType CenterOfRotationType;
  CenterOfRotationType m_CenterOfRotation;

  /** Test whether the specified tube point is inside the Image.
   * \param inputPoint the non-transformed tube point.
   * \param outputPoint the transformed tube point.
   * \param transform the transform to apply to the input point. */
  bool IsInside( const InputPointType & inputPoint,
    OutputPointType & outputPoint,
    const TransformType * transform ) const;

  ScalarType ComputeLaplacianMagnitude(
    const typename TubePointType::CovariantVectorType & tubeNormal,
    const ScalarType scale,
    const OutputPointType & currentPoint ) const;
  ScalarType ComputeThirdDerivatives(
    const CovariantVectorType & v,
    const ScalarType scale,
    const OutputPointType & currentPoint ) const;

  /**
   * \warning User is responsible for freeing the list, but not the elements
   * of the list.
   */
  typename TubeTreeType::ChildrenListType* GetTubes( void ) const;

  FeatureWeightsType m_FeatureWeights;
}; // End class ImageToTubeRigidMetric

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeImageToTubeRigidMetric.hxx"
#endif

#endif // End !defined( __itktubeImageToTubeRigidMetric_h )
