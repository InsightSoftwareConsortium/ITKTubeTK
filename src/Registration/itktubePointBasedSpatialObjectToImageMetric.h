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

#include <vector>

#include <itktubeSpatialObjectToImageMetric.h>

#include <itkTubeSpatialObject.h>
#include <itkSurfaceSpatialObject.h>
#include <itkPointBasedSpatialObject.h>

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

template< unsigned int ObjectDimension, class TFixedImage >
class PointBasedSpatialObjectToImageMetric
  : public SpatialObjectToImageMetric< ObjectDimension, TFixedImage >
{
public:
  /** Standard class typedefs. */
  typedef PointBasedSpatialObjectToImageMetric     Self;
  typedef SpatialObjectToImageMetric< ObjectDimension, TFixedImage >
                                                   Superclass;
  typedef SmartPointer< Self >                     Pointer;
  typedef SmartPointer< const Self >               ConstPointer;

  /**  Dimension of the image and tube.  */
  itkStaticConstMacro( ImageDimension, unsigned int,
    TFixedImage::ImageDimension );

  typedef TFixedImage                              FixedImageType;
  typedef SpatialObject< ObjectDimension >         MovingSpatialObjectType;

  /**  Type of the Transform Base class */
  typedef typename Superclass::TransformType       TransformType;
  typedef typename TransformType::Pointer          TransformPointer;
  typedef typename TransformType::ParametersType   TransformParametersType;
  typedef typename TransformType::JacobianType     TransformJacobianType;

  typedef typename Superclass::MovingPointType     MovingPointType;
  typedef typename Superclass::FixedPointType      FixedPointType;

  typedef PointBasedSpatialObject< ObjectDimension >
                                                   PointBasedSpatialObjectType;
  typedef typename PointBasedSpatialObjectType::SpatialObjectPointType
                                                   PointBasedSpatialObjectPointType;

  typedef PointBasedSpatialObject< ObjectDimension >
                                                   BlobType;
  typedef typename BlobType::SpatialObjectPointType
                                                   BlobPointType;
  typedef typename BlobType::SpatialObjectPointListType
                                                   BlobPointListType;
  typedef std::vector< double >                    BlobPointWeightListType;

  typedef TubeSpatialObject< ObjectDimension >     TubeType;
  typedef typename TubeType::TubePointType         TubePointType;
  typedef typename TubeType::TubePointListType     TubePointListType;
  typedef std::vector< double >                    TubePointWeightListType;

  typedef SurfaceSpatialObject< ObjectDimension >    SurfaceType;
  typedef typename SurfaceType::SurfacePointType     SurfacePointType;
  typedef typename SurfaceType::SurfacePointListType SurfacePointListType;
  typedef std::vector< double >                      SurfacePointWeightListType;

  typedef typename Superclass::ParametersType      ParametersType;
  typedef typename Superclass::MeasureType         MeasureType;
  typedef typename Superclass::FixedVectorType     FixedVectorType;
  typedef typename Superclass::MovingVectorType    MovingVectorType;
  typedef typename Superclass::DerivativeType      DerivativeType;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( PointBasedSpatialObjectToImageMetric,
    SpatialObjectToImageMetric );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Initialize the metric */
  void Initialize( void ) override;

  /** Get the Value for SingleValue Optimizers */
  MeasureType GetValue( const ParametersType & parameters ) const override;

  /** Get the Derivatives of the Match Measure */
  void GetDerivative( const ParametersType & parameters,
    DerivativeType & derivative ) const override;

  /** Get the Value And Derivative for SingleValue Optimizers */
  void GetValueAndDerivative( const ParametersType & parameters,
    MeasureType & Value, DerivativeType  & Derivative ) const override;

  /** Set/Get the portion of points used to compute the match metric. */
  itkSetMacro( SamplingRatio, double );
  itkGetConstMacro( SamplingRatio, double );

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

  void ComputeSubsampledPoints( void );
  void ComputeSubsampledPointsWeights( void );

  //itkSetMacro( SubsampledBlobPoints, BlobPointListType );
  itkGetConstMacro( SubsampledBlobPoints, BlobPointListType );
  //itkSetMacro( SubsampledBlobPointsWeights, BlobPointWeightListType );
  itkGetConstMacro( SubsampledBlobPointsWeights, BlobPointWeightListType );
  
  //itkSetMacro( SubsampledTubePoints, TubePointListType );
  itkGetConstMacro( SubsampledTubePoints, TubePointListType );
  //itkSetMacro( SubsampledTubePointsWeights, TubePointWeightListType );
  itkGetConstMacro( SubsampledTubePointsWeights, TubePointWeightListType );
  
  //itkSetMacro( SubsampledSurfacePoints, SurfacePointListType );
  itkGetConstMacro( SubsampledSurfacePoints, SurfacePointListType );
  //itkSetMacro( SubsampledSurfacePointsWeights, SurfacePointWeightListType );
  itkGetConstMacro( SubsampledSurfacePointsWeights, SurfacePointWeightListType );
  
  TransformPointer GetTransform( void ) const
    { return dynamic_cast<TransformType*>( this->m_Transform.GetPointer() ); }

  /** Get the number of transform parameters */
  inline unsigned int GetNumberOfParameters( void ) const override
    { return this->m_Transform->GetNumberOfParameters(); }

protected:

  // purposely not implemented
  PointBasedSpatialObjectToImageMetric( const Self& );
  void operator=( const Self& );

  unsigned int GetMaximumNumberOfPoints( void );

  virtual void ComputeCenterOfRotation( void );

  bool IsValidMovingPoint( const TubePointType & inputPoint ) const;
  bool IsValidMovingPoint( const SurfacePointType & inputPoint ) const;
  bool IsValidMovingPoint( const BlobPointType & inputPoint ) const;
  bool IsValidFixedPoint( const FixedPointType & fixedPoint ) const;

  typename MovingSpatialObjectType::ChildrenConstListType *
    GetPointBasedChildren( const MovingSpatialObjectType * parentSO,
      typename MovingSpatialObjectType::ChildrenConstListType * childrenSO=nullptr ) const;

private:

  PointBasedSpatialObjectToImageMetric( void );
  virtual ~PointBasedSpatialObjectToImageMetric( void );

  unsigned int m_MaximumNumberOfPoints;
  unsigned int m_MaximumNumberOfTubePoints;
  unsigned int m_MaximumNumberOfSurfacePoints;

  double     m_SamplingRatio;

  double     m_Kappa;
  double     m_Extent;

  double     m_TubePriorityRadius;
  double     m_TubeSamplingRadiusMin;
  double     m_TubeSamplingRadiusMax;

  MovingPointType  m_CenterOfRotation;

  /** Points with no tangets or normals */
  BlobPointListType          m_SubsampledBlobPoints;
  BlobPointWeightListType    m_SubsampledBlobPointsWeights;

  /** Points with one tangent and two normals */
  TubePointListType          m_SubsampledTubePoints;
  TubePointWeightListType    m_SubsampledTubePointsWeights;

  /** Points with one normal */
  SurfacePointListType       m_SubsampledSurfacePoints;
  SurfacePointWeightListType m_SubsampledSurfacePointsWeights;

}; // End class PointBasedSpatialObjectToImageMetric

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubePointBasedSpatialObjectToImageMetric.hxx"
#endif

#endif // End !defined( __itktubePointBasedSpatialObjectToImageMetric_h )
