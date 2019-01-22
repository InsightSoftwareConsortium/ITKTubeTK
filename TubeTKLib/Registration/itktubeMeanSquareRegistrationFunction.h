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

#ifndef __itktubeMeanSquareRegistrationFunction_h
#define __itktubeMeanSquareRegistrationFunction_h

#include <itkCentralDifferenceImageFunction.h>
#include <itkCovariantVector.h>
#include <itkInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkPDEDeformableRegistrationFunction.h>
#include <itkPoint.h>
#include <mutex>

namespace itk
{

namespace tube
{

/** \class MeanSquareRegistrationFunction
 *
 * This class encapsulate the PDE which drives the demons registration
 * algorithm. It is used by MeanSquareRegistrationFilter to compute the
 * output deformation field which will map a moving image onto a
 * a fixed image.
 *
 * Non-integer moving image values are obtained by using
 * interpolation. The default interpolator is of type
 * LinearInterpolateImageFunction. The user may set other
 * interpolators via method SetMovingImageInterpolator. Note that the input
 * interpolator must derive from base class InterpolateImageFunction.
 *
 * This class is templated over the fixed image type, moving image type,
 * and the deformation field type.
 *
 * \warning This filter assumes that the fixed image type, moving image type
 * and deformation field type all have the same number of dimensions.
 *
 * \sa MeanSquareRegistrationFilter
 * \ingroup FiniteDifferenceFunctions
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
class MeanSquareRegistrationFunction
  : public PDEDeformableRegistrationFunction< TFixedImage, TMovingImage,
                                              TDeformationField >
{
public:
  /** Standard class typedefs. */
  typedef MeanSquareRegistrationFunction    Self;
  typedef PDEDeformableRegistrationFunction< TFixedImage,
    TMovingImage, TDeformationField >       Superclass;
  typedef SmartPointer< Self >              Pointer;
  typedef SmartPointer< const Self >        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( MeanSquareRegistrationFunction,
    PDEDeformableRegistrationFunction );

  /** MovingImage image type. */
  typedef typename Superclass::MovingImageType     MovingImageType;
  typedef typename Superclass::MovingImagePointer  MovingImagePointer;
  typedef typename MovingImageType::PixelType      MovingImagePixelType;

  /** FixedImage image type. */
  typedef typename Superclass::FixedImageType     FixedImageType;
  typedef typename Superclass::FixedImagePointer  FixedImagePointer;
  typedef typename FixedImageType::IndexType      IndexType;
  typedef typename FixedImageType::SizeType       SizeType;
  typedef typename FixedImageType::SpacingType    SpacingType;

  /** Deformation field type. */
  typedef typename Superclass::DisplacementFieldType
    DeformationFieldType;
  typedef typename DeformationFieldType::Pointer
    DeformationFieldPointer;
  typedef typename DeformationFieldType::PixelType
    DeformationFieldPixelType;

  /** Inherit some enums from the superclass. */
  itkStaticConstMacro( ImageDimension, unsigned int,
    Superclass::ImageDimension );

  /** Inherit some enums from the superclass. */
  typedef typename Superclass::PixelType        PixelType;
  typedef typename Superclass::RadiusType       RadiusType;
  typedef typename Superclass::NeighborhoodType NeighborhoodType;
  typedef typename Superclass::FloatOffsetType  FloatOffsetType;
  typedef typename Superclass::TimeStepType     TimeStepType;

  /** Interpolator type. */
  typedef double                                     CoordRepType;
  typedef InterpolateImageFunction<MovingImageType, CoordRepType>
                                                     InterpolatorType;
  typedef typename InterpolatorType::Pointer         InterpolatorPointer;
  typedef typename InterpolatorType::PointType       PointType;
  typedef LinearInterpolateImageFunction<MovingImageType, CoordRepType>
                                                     DefaultInterpolatorType;

  /** Covariant vector type. */
  typedef CovariantVector< double, itkGetStaticConstMacro( ImageDimension ) >
    CovariantVectorType;

  /** Gradient calculator type. */
  typedef CentralDifferenceImageFunction<FixedImageType>
    GradientCalculatorType;
  typedef typename GradientCalculatorType::Pointer
    GradientCalculatorPointer;

  /** Set the moving image interpolator. */
  void SetMovingImageInterpolator( InterpolatorType * ptr )
    { m_MovingImageInterpolator = ptr; }

  /** Get the moving image interpolator. */
  InterpolatorType * GetMovingImageInterpolator( void )
    { return m_MovingImageInterpolator; }

  /** This class uses a constant time step of 1. */
  virtual TimeStepType ComputeGlobalTimeStep( void * itkNotUsed(
    globalData ) ) const
    { return m_TimeStep; }

  /** Return a pointer to a global data structure that is passed to
   * this object from the solver at each calculation.  */
  virtual void *GetGlobalDataPointer( void ) const
    {
    GlobalDataStruct *global = new GlobalDataStruct();
    return global;
    }

  /** Release memory for global data structure. */
  virtual void ReleaseGlobalDataPointer( void *GlobalData ) const
    { delete ( GlobalDataStruct * ) GlobalData;  }

  /** Set the object's state before each iteration. */
  virtual void InitializeIteration( void );

  /** This method is called by a finite difference solver image filter at
   * each pixel that does not lie on a data set boundary */
  virtual PixelType  ComputeUpdate( const NeighborhoodType &neighborhood,
    void *globalData,
    const FloatOffsetType &offset = FloatOffsetType( 0.0 ) );

  /** Computes the intensity difference between the fixed and moving image
   *  at the given index, under the given deformation vector. */
  virtual double ComputeIntensityDifference( const IndexType & index,
    const DeformationFieldPixelType & itvec );

  /** Get the energy mutex lock  */
  void SetEnergy( double energy )
    {
    m_EnergyCalculationLock.lock();
    this->m_Energy = energy;
    m_EnergyCalculationLock.unlock();
    }

  void SetBackgroundIntensity( MovingImagePixelType intensity )
    { m_BackgroundIntensity = intensity; }
  MovingImagePixelType GetBackgroundIntensity( void ) const
    { return m_BackgroundIntensity; }

  void SetIntensityDifferenceThreshold( double threshold )
    { m_IntensityDifferenceThreshold = threshold; }
  double GetIntensityDifferenceThreshold( void ) const
    { return m_IntensityDifferenceThreshold; }


protected:
  MeanSquareRegistrationFunction( void );
  ~MeanSquareRegistrationFunction( void ) {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

  /** FixedImage image neighborhood iterator type. */
  typedef ConstNeighborhoodIterator<FixedImageType>
    FixedImageNeighborhoodIteratorType;

  /** A global data type for this class of equation. Used to store
   * iterators for the fixed image. */
  struct GlobalDataStruct
    {
    FixedImageNeighborhoodIteratorType   m_FixedImageIterator;
    }; // End struct GlobalDataStruct

private:
  MeanSquareRegistrationFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  /** Cache fixed image information. */
  SpacingType                     m_FixedImageSpacing;

  /** Function to compute derivatives of the fixed image. */
  GradientCalculatorPointer       m_FixedImageGradientCalculator;

  /** Function to interpolate the moving image. */
  InterpolatorPointer             m_MovingImageInterpolator;

  /** The global time step. */
  TimeStepType                    m_TimeStep;

  /** Threshold below which the denominator term is considered zero. */
  double                          m_DenominatorThreshold;

  /** Threshold below which two intensity value are assumed to match. */
  double                          m_IntensityDifferenceThreshold;

  mutable std::mutex              m_EnergyCalculationLock;

  MovingImagePixelType            m_BackgroundIntensity;

}; // End class MeanSquareRegistrationFunction

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeMeanSquareRegistrationFunction.hxx"
#endif

#endif // End !defined( __itktubeMeanSquareRegistrationFunction_h )
