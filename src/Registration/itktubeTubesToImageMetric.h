/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#ifndef __itktubeTubesToImageMetric_h
#define __itktubeTubesToImageMetric_h

#include <itkCompensatedSummation.h>
#include <itkGaussianDerivativeImageFunction.h>
#include <itktubeSpatialObjectToImageMetric.h>

namespace itk
{

namespace tube
{

/** \class TubesToImageMetric
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

template< class TMovingTubeSpatialObject,
          class TFixedImage >
class TubesToImageMetric
  : public SpatialObjectToImageMetric< TMovingSpatialObject, TFixedImage >
{
public:
  /** Standard class typedefs. */
  typedef TubesToImageMetric                    Self;
  typedef SpatialObjectToImageMetric< TMovingSpatialObject, TFixedImage >
                                                Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;

  /**  Dimension of the image and tube.  */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TFixedImage::ImageDimension );
  itkStaticConstMacro( TubeDimension, unsigned int,
    TTubeSpatialObject::ObjectDimension );

  typedef TFixedImage                           FixedImageType;
  typedef TMovingSpatialObject                  GroupSpatialObjectType;
  typedef TTubeSpatialObject                    TubeType;
  typedef typename TubeType::TubePointType      TubePointType;

  typedef double                                ScalarType;
  typedef GaussianDerivativeImageFunction< TFixedImage >
                                                DerivativeImageFunctionType;
  typedef typename Superclass::DerivativeType   DerivativeType;
  typedef typename Superclass::ParametersType   ParametersType;
  typedef typename Superclass::MeasureType      MeasureType;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( TubesToImageMetric, ImageToSpatialObjectMetric );

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
  typedef TTransform                               TransformType;
  typedef typename TransformType::Pointer          TransformPointer;
  typedef typename TransformType::InputPointType   InputPointType;
  typedef typename TransformType::OutputPointType  OutputPointType;
  typedef typename TransformType::ParametersType   TransformParametersType;
  typedef typename TransformType::JacobianType     TransformJacobianType;

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

  /** The radius mulitplier used to define the scale of the Gaussian kernel
   *  used to make image measures in the metric. */
  itkSetMacro( Kappa, ScalarType );
  itkGetConstMacro( Kappa, ScalarType );

  /** Set/Get the extent of the Gaussian kernel used for image measures. */
  itkSetMacro( Extent, ScalarType );
  itkGetConstMacro( Extent, ScalarType );

  /** Set/Get the minimum tube radius used by the metric. */
  itkSetMacro( SamplingRadiusMin, ScalarType );
  itkGetConstMacro( SamplingRadiusMin, ScalarType );

  /** Set/Get the maximum tube radius used by the metric. */
  itkSetMacro( SamplingRadiusMax, ScalarType );
  itkGetConstMacro( SamplingRadiusMax, ScalarType );

  /** Set/Get the space between points at which the metric is computed. Spacing
   * is given as a multiplier of the radius of the adjacent points. */
  itkSetMacro( SampleSpacingInRadii, ScalarType );
  itkGetConstMacro( SampleSpacingInRadii, ScalarType );

  void ResampleTube( void );
  
  TransformPointer GetTransform( void ) const
    {
    return dynamic_cast<TransformType*>( this->m_Transform.GetPointer() );
    }

  /** Downsample the tube points by this integer value. */

protected:
  TubesToImageMetric( void );
  virtual ~TubesToImageMetric( void );

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
  TubesToImageMetric( const Self& ); // purposely not implemented
  void operator=( const Self& ); // purposely not implemented

  typename DerivativeImageFunctionType::Pointer m_DerivativeImageFunction;

  ScalarType m_Kappa;
  ScalarType m_Extent;
  ScalarType m_SamplingRadiusMin;
  ScalarType m_SamplingRadiusMax;
  ScalarType m_SamplingSpacingInRadii;

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

}; // End class TubesToImageMetric

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTubesToImageMetric.hxx"
#endif

#endif // End !defined( __itktubeTubesToImageMetric_h )
