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
#ifndef __itkImageToTubeRigidMetric_h
#define __itkImageToTubeRigidMetric_h

#include "itkPoint.h"
#include "itkGroupSpatialObject.h"
#include "itkVesselTubeSpatialObject.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkEuler3DTransform.h"
#include "itkImageToSpatialObjectMetric.h"
#include <itkGaussianDerivativeImageFunction.h>

namespace itk
{

/**
 * \class ImageToTubeRigidMetric
 * \brief Computes similarity between two objects to be registered
 * The metric implemented here corresponds to the following paper:
 * \link http://www.cs.unc.edu/Research/MIDAG/pubs/papers/MICCAI01-aylwardVReg.pdf
 * The metric is based on the fact that vessel centerlines are scaled
 * intensity ridges in the image.
 * \warning (Derivative)
*/

template < class TFixedImage, class TMovingTube>
class ITK_EXPORT ImageToTubeRigidMetric
: public ImageToSpatialObjectMetric<TFixedImage, TMovingTube>
{
public:
  typedef ImageToTubeRigidMetric                Self;
  typedef ImageToSpatialObjectMetric<TFixedImage, TMovingTube>
                                                Superclass;
  typedef SmartPointer<Self>                    Pointer;
  typedef SmartPointer<const Self>              ConstPointer;

  /**  Dimension of the image and tube.  */
  itkStaticConstMacro( ImageDimension, unsigned int, TFixedImage::ImageDimension );
  itkStaticConstMacro( TubeDimension, unsigned int, TMovingTube::ObjectDimension );

  /** Type definition for a tube point */
  typedef VesselTubeSpatialObjectPoint< TubeDimension>  TubePointType;
  typedef VesselTubeSpatialObject< TubeDimension >      TubeType;
  typedef TMovingTube                                   MovingTubeType;
  typedef typename MovingTubeType::ChildrenListType
                                                        ChildrenListType;
  typedef GroupSpatialObject< TubeDimension >           TubeNetType;
  typedef TFixedImage                                   FixedImageType;
  typedef typename FixedImageType::IndexType            IndexType;
  typedef GaussianDerivativeImageFunction< TFixedImage >
                                                DerivativeImageFunctionType;

  typedef double                                InternalComputationValueType;
  typedef typename Superclass::DerivativeType   DerivativeType;
  typedef typename Superclass::ParametersType   ParametersType;
  typedef typename Superclass::MeasureType      MeasureType;

  typedef vnl_vector< InternalComputationValueType >                    VectorType;
  typedef vnl_matrix< InternalComputationValueType >                    MatrixType;
  typedef Point< InternalComputationValueType, ImageDimension >         PointType;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ImageToTubeRigidMetric, ImageToSpatialObjectMetric );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Space dimension is the dimension of parameters space */
  enum { SpaceDimension = 6 };

  enum { RangeDimension = 6 };

  unsigned int GetNumberOfParameters( void ) const
    { return SpaceDimension; }

  /** Typedef for the Range calculator */
  typedef MinimumMaximumImageCalculator<FixedImageType> RangeCalculatorType;

  /** Type used for representing point components  */
  typedef typename Superclass::CoordinateRepresentationType
                                               CoordinateRepresentationType;

  /** Type definition for the size */
  typedef typename TFixedImage::SizeType       SizeType;

  /** Type definition for the pixel type */
  typedef typename TFixedImage::PixelType      PixelType;

  /**  Type of the Transform Base class */
  typedef Euler3DTransform<double>                 TransformType;
  typedef typename TransformType::Pointer          TransformPointer;
  typedef typename TransformType::InputPointType   InputPointType;
  typedef typename TransformType::OutputPointType  OutputPointType;
  typedef typename TransformType::ParametersType   TransformParametersType;
  typedef typename TransformType::JacobianType     TransformJacobianType;

  /** Get the Derivatives of the Match Measure */
  const DerivativeType & GetDerivative( const ParametersType &
    parameters ) const;
  void GetDerivative( const ParametersType & parameters,
    DerivativeType & derivative ) const;

  /** Get the Value for SingleValue Optimizers */
  MeasureType  GetValue( const ParametersType & parameters ) const;

  /** Get Value and Derivatives for MultipleValuedOptimizers */
  void GetValueAndDerivative( const ParametersType & parameters,
    MeasureType & Value, DerivativeType  & Derivative ) const;

  /** SubSample the MovingSpatialObject tube */
  void SubSampleTube();

  /** Apply the center of rotation to the transformation */
  ParametersType ApplyCenterOfRotation( const ParametersType & parameters );

  /** Set kappa value */
  itkSetMacro( Kappa, double );

  vnl_vector_fixed<double, 3> GetCenterOfRotation( void )
    { return m_RotationCenter; }

  /** Initialize the metric */
  void Initialize( void ) throw ( ExceptionObject );

  /** Set the extent of the blurring */
  itkSetMacro( Extent, double );
  itkGetMacro( Extent, double );

  TransformPointer GetTransform( void ) const
    { return dynamic_cast<TransformType*>( this->m_Transform.GetPointer() ); }

  /** Downsample the tube points by this integer value. */
  itkSetMacro( Sampling, OffsetValueType );
  itkGetMacro( Sampling, OffsetValueType );

protected:
  ImageToTubeRigidMetric();
  virtual ~ImageToTubeRigidMetric();

  ImageToTubeRigidMetric( const Self& ); // purposely not implemented
  void operator=( const Self& ); // purposely not implemented

  void ComputeImageRange( void );

  void GetDeltaAngles( const Point< InternalComputationValueType, 3> & x,
                       const vnl_vector_fixed< InternalComputationValueType, 3> & dx,
                       const vnl_vector_fixed< InternalComputationValueType, 3> & offsets,
                       double angle[3] ) const;

private:
  typename DerivativeImageFunctionType::Pointer m_DerivativeImageFunction;

  unsigned int                               m_NumberOfPoints;
  std::list< InternalComputationValueType >  m_Weight;
  InternalComputationValueType               m_SumWeight;
  InternalComputationValueType               m_ImageMin;
  InternalComputationValueType               m_ImageMax;
  typename RangeCalculatorType::Pointer      m_RangeCalculator;
  unsigned int                               m_Iteration;
  double                                     m_Kappa;
  double                                     m_Extent;
  InternalComputationValueType               m_InitialScale;
  OffsetValueType                            m_Sampling;

  mutable OutputPointType                m_CurrentPoint;
  mutable InternalComputationValueType   m_BlurredValue;

  vnl_vector_fixed< InternalComputationValueType, TubeDimension >  m_Offsets;
  vnl_vector_fixed< InternalComputationValueType, 3 >              m_RotationCenter;

  vnl_vector_fixed< InternalComputationValueType, 3 >              m_Factors;

  // TODO Create gfact var insteal calculating it all the time

  /** Set the scale of the blurring */
  itkGetConstMacro( InitialScale, InternalComputationValueType );

  /** Test whether the specified point is inside
  This method overload the one in the ImageMapper class
  \warning This method cannot be safely used in more than one thread at a time.
  \sa Evaluate() */
  bool IsInside( const InputPointType & point ) const;

  VectorType *  EvaluateAllDerivatives( void ) const;

  itkGetConstMacro( BlurredValue, InternalComputationValueType );

  double ComputeLaplacianMagnitude( Vector< InternalComputationValueType, 3 > *v,
    const InternalComputationValueType & scale ) const;
  double ComputeThirdDerivatives( Vector< InternalComputationValueType, 3 > *v,
    const InternalComputationValueType & scale ) const;
  double ComputeDerivatives( Vector< InternalComputationValueType, 3 > *v ) const;

  void ClampPointBoundsToImage( int bounds[6] );
  void ComputeCenterRotation();
  typename TubeNetType::ChildrenListType* GetTubes() const;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToTubeRigidMetric.txx"
#endif

#endif
