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
#ifndef __itkImageToTubeRigidMetric2_h
#define __itkImageToTubeRigidMetric2_h

#include "itkPoint.h"
#include "itkGroupSpatialObject.h"
#include "itkVesselTubeSpatialObject.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkLinearInterpolateImageFunction.h"
#include <time.h>
#include "itkEuler3DTransform.h"
#include "itkImageToSpatialObjectMetric.h"
#include <vnl/vnl_vector.h>
#include <itkGaussianDerivativeImageFunction.h>

namespace itk
{
/**
 * \class ImageToTubeRigidMetric2
 * \brief Computes similarity between two objects to be registered
 * The metric implemented here corresponds to the following paper:
 * \link http://www.cs.unc.edu/Research/MIDAG/pubs/papers/MICCAI01-aylwardVReg.pdf
 * The metric is based on the fact that vessel centerlines are scaled
 * intensity ridges in the image. The improvment implemented here is about
 * the precomputation of the internal kernels proposed by Lange et al.:
 * \link http://www.zib.de/lamecker/publications/cars2007.pdf
*/

template < class TFixedImage, class TMovingSpatialObject>
class ITK_EXPORT ImageToTubeRigidMetric2
: public ImageToSpatialObjectMetric<TFixedImage, TMovingSpatialObject>
{
public:
  typedef ImageToTubeRigidMetric2                 Self;
  typedef ImageToSpatialObjectMetric<TFixedImage, TMovingSpatialObject>
                                                  Superclass;
  typedef SmartPointer<Self>                      Pointer;
  typedef SmartPointer<const Self>                ConstPointer;
  typedef VesselTubeSpatialObjectPoint<3>         TubePointType;
  typedef VesselTubeSpatialObject<3>              TubeType;
  typedef TMovingSpatialObject                    MovingSpatialObjectType;
  typedef typename MovingSpatialObjectType::ChildrenListType
                                                  ChildrenListType;
  typedef GroupSpatialObject<3>                   TubeNetType;
  typedef Image<unsigned char, 3>                 MaskImageType;
  typedef typename MaskImageType::Pointer         MaskImagePointer;
  typedef typename MaskImageType::IndexType       IndexType;
  typedef TFixedImage                             FixedImageType;
  typedef GaussianDerivativeImageFunction<TFixedImage>
                                                  DerivativeImageFunctionType;
  typedef typename Superclass::DerivativeType     DerivativeType;
  typedef typename Superclass::ParametersType     ParametersType;
  typedef typename Superclass::MeasureType        MeasureType;
  typedef vnl_vector<double>                      VectorType;
  typedef vnl_matrix<double>                      MatrixType;
  typedef Point<double, 3>                        PointType;
  typedef std::list<PointType>                    PointsType;

  struct TransformedPointType
    {
    TransformedPointType(double s = 1.0, double w = 0.0) :
      Scale( s ), Weight( w ) {}

    TubePointType TransformedPoint;
    double Scale;
    double Weight;
    };
  typedef std::vector<TransformedPointType> TransformedPointsType;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ImageToTubeRigidMetric2, ImageToSpatialObjectMetric );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Space dimension is the dimension of parameters space */
  enum { SpaceDimension = 6 };
  enum { ImageDimension = 3 };
  enum { RangeDimension = 6 };

  unsigned int GetNumberOfParameters( void ) const { return SpaceDimension; }

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
  const DerivativeType& GetDerivative( const ParametersType&
                                        parameters ) const;
  void GetDerivative( const ParametersType & parameters,
                      DerivativeType& derivative ) const;

  /** Get the Value for SingleValue Optimizers */
  MeasureType  GetValue( const ParametersType & parameters ) const;

  /** Get Value and Derivatives for MultipleValuedOptimizers */
  void GetValueAndDerivative( const ParametersType& parameters,
                              MeasureType& Value,
                              DerivativeType& Derivative ) const;

  /** Apply the center of rotation to the transformation */
  ParametersType ApplyCenterOfRotation( const ParametersType & parameters );

  /** Set kappa value, control the radius scalling  of the metric */
  itkSetMacro( Kappa, double );

  /** Return the center of rotation of the whole extent */
  vnl_vector_fixed<double, 3> GetCenterOfRotation( void )
    { return m_RotationCenter; }

  /** Initialize the metric */
  void Initialize( void ) throw ( ExceptionObject );

  /** Set the extent of the blurring */
  itkSetMacro( Extent, double );
  itkGetMacro( Extent, double );

  /** Return the current transform */
  TransformPointer GetTransform( void ) const
    { return dynamic_cast<TransformType*>( this->m_Transform.GetPointer() ); }

  /** Set the verbosity of the computations */
  itkSetMacro( Verbose, bool );
  itkGetMacro( Verbose, bool );

  /** Control the ratio of the subsampling */
  itkSetMacro( Sampling, unsigned int );
  itkGetMacro( Sampling, unsigned int );

protected:

  ImageToTubeRigidMetric2();
  virtual ~ImageToTubeRigidMetric2(){}
  ImageToTubeRigidMetric2( const Self& ){}
  void operator=( const Self& ){}

  void ComputeImageRange( void );
  void GetDeltaAngles( const Point<double, 3> & x,
                       const vnl_vector_fixed<double, 3> & dx,
                       double angle[3] ) const;

private:

  typename DerivativeImageFunctionType::Pointer m_DerivativeImageFunction;

  unsigned int                           m_NumberOfPoints;
  mutable double                         m_SumWeight;
  double                                 m_ImageMin;
  double                                 m_ImageMax;
  typename RangeCalculatorType::Pointer  m_RangeCalculator;
  unsigned int                           m_Iteration;
  double                                 m_Kappa;
  double                                 m_Extent;
  unsigned int                           m_Sampling;
  mutable double                         m_CachedValue;
  mutable DerivativeType                 m_CachedDerivative;
  bool                                   m_Verbose;
  mutable vnl_vector_fixed<double, 3>    m_Offset;
  vnl_vector_fixed<double, 3>            m_RotationCenter;
  vnl_vector_fixed<double, 3>            m_Factors;
  mutable TransformedPointsType          m_TransformedPoints;

  void SubSampleTube();
  void SetOffset( double oX, double oY, double oZ ) const;
  void ComputeCenterRotation();
  void UpdateTransformedPoints() const;
  bool UpdateTransformedPoints( const ParametersType & parameters ) const;

  TubeNetType::ChildrenListType* GetTubes() const;

  double ComputeLaplacianMagnitude( TransformedPointType& ) const;
  double ComputeThirdDerivatives( Vector<double, 3> *v,
                                  TransformedPointType& ) const;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToTubeRigidMetric2.txx"
#endif

#endif
