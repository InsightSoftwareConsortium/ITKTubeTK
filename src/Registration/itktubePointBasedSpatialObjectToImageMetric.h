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

#ifndef __itktubePointBasedSpatialObjectToImageMetric_h
#define __itktubePointBasedSpatialObjectToImageMetric_h

#include <itkCompensatedSummation.h>
#include <itktubeSpatialObjectToImageMetric.h>

namespace itk
{

namespace tube
{

/** \class PointBasedSpatialObjectToImageMetric
 * \brief Computes similarity between two objects to be registered
 * The metric implemented here corresponds to the following paper:
 * \link
 * http://www.cs.unc.edu/Research/MIDAG/pubs/papers/MICCAI01-aylwardVReg.pdf
 * The metric is based on the fact that vessel centerlines are scaled
 * intensity ridges in the image.
 *
 * \warning ( Derivative )
 */

template< unsigned int ObjectDimension, class TFixedImage, class  >
class PointBasedSpatialObjectToImageMetric
  : public SpatialObjectToImageMetric< ObjectDimension, TFixedImage >
{
public:
  /** Standard class typedefs. */
  typedef PointBasedSpatialObjectToImageMetric  Self;
  typedef SpatialObjectToImageMetric< ObjectDimension, TFixedImage >
                                                Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;

  /**  Dimension of the image and tube.  */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TFixedImage::ImageDimension );

  typedef TFixedImage                           FixedImageType;
  typedef SpatialObject< ObjectDimension >      SpatialObjectType;

  typedef PointBasedSpatialObject< ObjectDimension >
    PointBasedSpatialObjectType;
  typedef typename PointBasedSpatialObjectType::PointType
    SpatialObjectPointType;
  typedef typename PointBasedSpatialObjectType::SpatialObjectPointListType
    SpatialObjectPointListType;

  typedef TubeSpatialObject< ObjectDimension >  TubeType;
  typedef typename TubeType::PointType          TubePointType;
  typedef typename TubeType::TubePointListType  TubePointListType;

  typedef SurfacSpatialObject< ObjectDimension >     SurfaceType;
  typedef typename SurfaceType::PointType            SurfacePointType;
  typedef typename SurfaceType::SurfacePointListType SurfacePointListType;

  typedef typename Superclass::DerivativeType   DerivativeType;
  typedef typename Superclass::ParametersType   ParametersType;
  typedef typename Superclass::MeasureType      MeasureType;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( PointBasedSpatialObjectToImageMetric,
    SpatialObjectToImageMetric );

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
  typedef typename Superclass::TransformType       TransformType;
  typedef typename TransformType::Pointer          TransformPointer;
  typedef typename TransformType::InputPointType   InputPointType;
  typedef typename TransformType::OutputPointType  OutputPointType;
  typedef typename TransformType::ParametersType   TransformParametersType;
  typedef typename TransformType::JacobianType     TransformJacobianType;

  /** Initialize the metric */
  void Initialize( void ) override;

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

  /** The radius (for tube points) or image voxel size (for all other points)
   *    mulitplier used to define the scale of the Gaussian kernel
   *    used to make image measures in the metric. */
  itkSetMacro( Kappa, double );
  itkGetConstMacro( Kappa, double );

  /** Set/Get the extent of the Gaussian kernel used for image measures. */
  itkSetMacro( Extent, double );
  itkGetConstMacro( Extent, double );

  /** Set/Get the minimum tube radius used by the metric. */
  itkSetMacro( TubeSamplingRadiusMin, double );
  itkGetConstMacro( TubeSamplingRadiusMin, double );

  /** Set/Get the maximum tube radius used by the metric. */
  itkSetMacro( TubeSamplingRadiusMax, double );
  itkGetConstMacro( TubeSamplingRadiusMax, double );

  /** Set/Get the priority tube radius used in computing weights. */
  itkSetMacro( TubePriorityRadius, double );
  itkGetConstMacro( TubePriorityRadius, double );

  void SubsamplePoints( void );
  void ComputeSubsampledPointsWeights( void );

  itkSetMacro( SubsampledPoints, SpatialObjectPointListType );
  itkGetConstMacro( SubsampledPoints, SpatialObjectPointListType );
  
  itkSetMacro( SubsampledTubePoints, TubePointListType );
  itkGetConstMacro( SubsampledTubePoints, TubePointListType );
  
  itkSetMacro( SubsampledSurfacePoints, SurfacePointListType );
  itkGetConstMacro( SubsampledSurfacePoints, SurfacePointListType );
  
  TransformPointer GetTransform( void ) const
    {
    return dynamic_cast<TransformType*>( this->m_Transform.GetPointer() );
    }

protected:
  PointBasedSpatialObjectToImageMetric( void );
  virtual ~PointBasedSpatialObjectToImageMetric( void );

  typedef CovariantVector< double, ObjectDimension >    CovariantVectorType;
  typedef Vector< double, ObjectDimension >             VectorType;

  typedef Matrix< double, ObjectDimension, ObjectDimension >   MatrixType; 

  typedef typename SpatialObjectPointType::PointType     PointType;

  typedef vnl_vector< double >                       VnlVectorType;
  typedef vnl_matrix< double >                       VnlMatrixType;

  unsigned int GetMaximumNumberOfPoints( void );

  virtual void ComputeCenterOfRotation( void );

private:
  // purposely not implemented
  PointBasedSpatialObjectToImageMetric( const Self& );
  void operator=( const Self& );

  double m_Kappa;
  double m_Extent;

  double m_TubePriorityRadius;
  double m_TubeSamplingRadiusMin;
  double m_TubeSamplingRadiusMax;

  /** The center of rotation of the weighted tube points. */
  PointType m_CenterOfRotation;

  /** Points with no tangets or normals */
  SpatialObjectPointListType m_SubsampledPoints;

  /** Points with one tangent and two normals */
  TubePointListType          m_SubsampledTubePoints;

  /** Points with one normal */
  SurfacePointListType       m_SubsampledSurfacePoints;

  /** Test whether the specified tube point is inside the Image.
   * \param inputPoint the non-transformed tube point.
   * \param outputPoint the transformed tube point.
   * \param transform the transform to apply to the input point. */
  bool IsValidMovingPoint( const TubePointType & inputPoint ) const;
  bool IsValidFixedPoint( OutputPointType & outputPoint ) const;

  /**
   * \warning User is responsible for freeing the list, but not the elements
   * of the list.
   */
  typename SpatialObjectType::ChildrenListType *
    GetPointBasedChildren( SpatialObjectType::Pointer & parentSO,
      SpatialObjectType::ChildrenListType * childrenSO=nullptr ) const;

}; // End class PointBasedSpatialObjectToImageMetric

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubePointBasedSpatialObjectToImageMetric.hxx"
#endif

#endif // End !defined( __itktubePointBasedSpatialObjectToImageMetric_h )
