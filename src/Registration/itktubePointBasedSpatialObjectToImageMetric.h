/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef itktubePointBasedSpatialObjectToImageMetric_h
#define itktubePointBasedSpatialObjectToImageMetric_h

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
 * https://www.cs.unc.edu/Research/MIDAG/pubs/papers/MICCAI01-aylwardVReg.pdf
 * The metric is based on the fact that vessel centerlines are scaled
 * intensity ridges in the image.
 *
 * \warning ( Derivative )
 */

template <unsigned int ObjectDimension, class TFixedImage>
class PointBasedSpatialObjectToImageMetric : public SpatialObjectToImageMetric<ObjectDimension, TFixedImage>
{
public:
  /** Standard class type alias. */
  using Self = PointBasedSpatialObjectToImageMetric;
  using Superclass = SpatialObjectToImageMetric<ObjectDimension, TFixedImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /**  Dimension of the image and tube.  */
  itkStaticConstMacro(ImageDimension, unsigned int, TFixedImage::ImageDimension);

  using FixedImageType = TFixedImage;
  using MovingSpatialObjectType = SpatialObject<ObjectDimension>;

  /**  Type of the Transform Base class */
  using TransformType = typename Superclass::TransformType;
  using TransformPointer = typename TransformType::Pointer;
  using TransformParametersType = typename TransformType::ParametersType;
  using TransformJacobianType = typename TransformType::JacobianType;

  using MovingPointType = typename Superclass::MovingPointType;
  using FixedPointType = typename Superclass::FixedPointType;

  using PointBasedSpatialObjectType = PointBasedSpatialObject<ObjectDimension>;
  typedef typename PointBasedSpatialObjectType::SpatialObjectPointType PointBasedSpatialObjectPointType;

  using BlobType = PointBasedSpatialObject<ObjectDimension>;
  typedef typename BlobType::SpatialObjectPointType     BlobPointType;
  typedef typename BlobType::SpatialObjectPointListType BlobPointListType;
  using BlobPointWeightListType = std::vector<double>;

  using TubeType = TubeSpatialObject<ObjectDimension>;
  using TubePointType = typename TubeType::TubePointType;
  using TubePointListType = typename TubeType::TubePointListType;
  using TubePointWeightListType = std::vector<double>;

  using SurfaceType = SurfaceSpatialObject<ObjectDimension>;
  using SurfacePointType = typename SurfaceType::SurfacePointType;
  using SurfacePointListType = typename SurfaceType::SurfacePointListType;
  using SurfacePointWeightListType = std::vector<double>;

  using ParametersType = typename Superclass::ParametersType;
  using MeasureType = typename Superclass::MeasureType;
  using FixedVectorType = typename Superclass::FixedVectorType;
  using MovingVectorType = typename Superclass::MovingVectorType;
  using DerivativeType = typename Superclass::DerivativeType;

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro(PointBasedSpatialObjectToImageMetric);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Initialize the metric */
  void
  Initialize(void) override;

  /** Get the Value for SingleValue Optimizers */
  MeasureType
  GetValue(const ParametersType & parameters) const override;

  /** Get the Derivatives of the Match Measure */
  void
  GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const override;

  /** Get the Value And Derivative for SingleValue Optimizers */
  void
  GetValueAndDerivative(const ParametersType & parameters,
                        MeasureType &          Value,
                        DerivativeType &       Derivative) const override;

  /** Set/Get the portion of points used to compute the match metric. */
  itkSetMacro(SamplingRatio, double);
  itkGetConstMacro(SamplingRatio, double);

  /** The radius (for tube points) or image voxel size (for all other points)
   *    mulitplier used to define the scale of the Gaussian kernel
   *    used to make image measures in the metric. */
  itkSetMacro(Kappa, double);
  itkGetConstMacro(Kappa, double);

  /** Set/Get the extent of the Gaussian kernel used for image measures. */
  itkSetMacro(Extent, double);
  itkGetConstMacro(Extent, double);

  /** Set/Get the minimum tube radius used by the metric. */
  itkSetMacro(TubeSamplingRadiusMin, double);
  itkGetConstMacro(TubeSamplingRadiusMin, double);
  /** Set/Get the maximum tube radius used by the metric. */
  itkSetMacro(TubeSamplingRadiusMax, double);
  itkGetConstMacro(TubeSamplingRadiusMax, double);
  /** Set/Get the priority tube radius used in computing weights. */
  itkSetMacro(TubePriorityRadius, double);
  itkGetConstMacro(TubePriorityRadius, double);

  void
  ComputeSubsampledPoints(void);
  void
  ComputeSubsampledPointsWeights(void);

  // itkSetMacro( SubsampledBlobPoints, BlobPointListType );
  itkGetConstMacro(SubsampledBlobPoints, BlobPointListType);
  // itkSetMacro( SubsampledBlobPointsWeights, BlobPointWeightListType );
  itkGetConstMacro(SubsampledBlobPointsWeights, BlobPointWeightListType);

  // itkSetMacro( SubsampledTubePoints, TubePointListType );
  itkGetConstMacro(SubsampledTubePoints, TubePointListType);
  // itkSetMacro( SubsampledTubePointsWeights, TubePointWeightListType );
  itkGetConstMacro(SubsampledTubePointsWeights, TubePointWeightListType);

  // itkSetMacro( SubsampledSurfacePoints, SurfacePointListType );
  itkGetConstMacro(SubsampledSurfacePoints, SurfacePointListType);
  // itkSetMacro( SubsampledSurfacePointsWeights, SurfacePointWeightListType );
  itkGetConstMacro(SubsampledSurfacePointsWeights, SurfacePointWeightListType);

  TransformPointer
  GetTransform(void) const
  {
    return dynamic_cast<TransformType *>(this->m_Transform.GetPointer());
  }

  /** Get the number of transform parameters */
  inline unsigned int
  GetNumberOfParameters(void) const override
  {
    return this->m_Transform->GetNumberOfParameters();
  }

protected:
  // purposely not implemented
  PointBasedSpatialObjectToImageMetric(const Self &);
  void
  operator=(const Self &);

  unsigned int
  GetMaximumNumberOfPoints(void);

  virtual void
  ComputeCenterOfRotation(void);

  bool
  IsValidMovingPoint(const TubePointType & inputPoint) const;
  bool
  IsValidMovingPoint(const SurfacePointType & inputPoint) const;
  bool
  IsValidMovingPoint(const BlobPointType & inputPoint) const;
  bool
  IsValidFixedPoint(const FixedPointType & fixedPoint) const;

  typename MovingSpatialObjectType::ChildrenConstListType *
  GetPointBasedChildren(const MovingSpatialObjectType *                           parentSO,
                        typename MovingSpatialObjectType::ChildrenConstListType * childrenSO = nullptr) const;

private:
  PointBasedSpatialObjectToImageMetric(void);
  virtual ~PointBasedSpatialObjectToImageMetric(void);

  unsigned int m_MaximumNumberOfPoints;
  unsigned int m_MaximumNumberOfTubePoints;
  unsigned int m_MaximumNumberOfSurfacePoints;

  double m_SamplingRatio;

  double m_Kappa;
  double m_Extent;

  double m_TubePriorityRadius;
  double m_TubeSamplingRadiusMin;
  double m_TubeSamplingRadiusMax;

  MovingPointType m_CenterOfRotation;

  /** Points with no tangets or normals */
  BlobPointListType       m_SubsampledBlobPoints;
  BlobPointWeightListType m_SubsampledBlobPointsWeights;

  /** Points with one tangent and two normals */
  TubePointListType       m_SubsampledTubePoints;
  TubePointWeightListType m_SubsampledTubePointsWeights;

  /** Points with one normal */
  SurfacePointListType       m_SubsampledSurfacePoints;
  SurfacePointWeightListType m_SubsampledSurfacePointsWeights;

}; // End class PointBasedSpatialObjectToImageMetric

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubePointBasedSpatialObjectToImageMetric.hxx"
#endif

#endif // End !defined( __itktubePointBasedSpatialObjectToImageMetric_h )
