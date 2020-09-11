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
#ifndef __tubeRegisterImages_h
#define __tubeRegisterImages_h

// ITK Includes
#include "itkProcessObject.h"

// TubeTK Includes
#include "tubeWrappingMacros.h"

#include "itkImageToImageRegistrationHelper.h"

namespace tube
{
/** \class RegisterImages
 *
 *  \ingroup TubeTK
 */

template< class TImage >
class RegisterImages:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef RegisterImages                  Self;
  typedef itk::ProcessObject              Superclass;
  typedef itk::SmartPointer< Self >       Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  typedef TImage                                             ImageType;

  typedef itk::ImageToImageRegistrationHelper< ImageType >   FilterType;

  typedef typename FilterType::MaskObjectType       MaskObjectType;
  typedef typename FilterType::PointType            PointType;
  typedef typename FilterType::PixelType            PixelType;

  typedef typename FilterType::LandmarkVectorType   LandmarkVectorType;

  typedef typename FilterType::RigidTransformType   RigidTransformType;
  typedef typename FilterType::AffineTransformType  AffineTransformType;

  typedef typename FilterType::MatrixTransformType  MatrixTransformType;
  typedef typename FilterType::BSplineTransformType BSplineTransformType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( RegisterImages, ProcessObject );

  /* Set input image */
  tubeWrapCallWithConstReferenceArgMacro( LoadFixedImage, std::string, Filter );
  tubeWrapSetConstObjectMacro( FixedImage, ImageType, Filter );

  tubeWrapCallWithConstReferenceArgMacro( LoadMovingImage, std::string, Filter );
  tubeWrapSetConstObjectMacro( MovingImage, ImageType, Filter );

  tubeWrapSetMacro( RandomNumberSeed, unsigned int, Filter );
  tubeWrapGetMacro( RandomNumberSeed, unsigned int, Filter );

  tubeWrapSetMacro( UseFixedImageMaskObject, bool, Filter );
  tubeWrapGetMacro( UseFixedImageMaskObject, bool, Filter );
  tubeWrapSetConstObjectMacro( FixedImageMaskObject, MaskObjectType, Filter );
  tubeWrapGetConstObjectMacro( FixedImageMaskObject, MaskObjectType, Filter );

  tubeWrapSetMacro( UseMovingImageMaskObject, bool, Filter );
  tubeWrapGetMacro( UseMovingImageMaskObject, bool, Filter );
  tubeWrapSetConstObjectMacro( MovingImageMaskObject, MaskObjectType, Filter );
  tubeWrapGetConstObjectMacro( MovingImageMaskObject, MaskObjectType, Filter );

  tubeWrapSetMacro( UseRegionOfInterest, bool, Filter );
  tubeWrapGetMacro( UseRegionOfInterest, bool, Filter );
  tubeWrapSetMacro( RegionOfInterestPoint1, PointType, Filter );
  tubeWrapGetMacro( RegionOfInterestPoint1, PointType, Filter );
  tubeWrapSetMacro( RegionOfInterestPoint2, PointType, Filter );
  tubeWrapGetMacro( RegionOfInterestPoint2, PointType, Filter );

  tubeWrapSetMacro( SampleFromOverlap, bool, Filter );
  tubeWrapGetMacro( SampleFromOverlap, bool, Filter );

  tubeWrapSetMacro( SampleIntensityPortion, double, Filter );
  tubeWrapGetMacro( SampleIntensityPortion, double, Filter );

  //
  // Update
  //
  tubeWrapUpdateMacro( Filter );

  tubeWrapCallMacro( Initialize, Filter );

  //
  // Resample
  //
  const ImageType * ResampleImage(
    std::string interp="LINEAR",
    const ImageType * movingImage = nullptr,
    const MatrixTransformType * matrixTransform = nullptr,
    //const BSplineTransformType * bsplineTransform = nullptr,
    PixelType defaultPixelValue = 0,
    double portionOfTransformToApply = 1.0 );

  const ImageType * GetFinalMovingImage(
    std::string interp="LINEAR_INTERPOLATION",
    PixelType defaultPixelValue = 0 );

  //
  // Compare results with a baseline
  //
  tubeWrapCallWithConstReferenceArgMacro( LoadBaselineImage,
    std::string, Filter );
  tubeWrapSetConstObjectMacro( BaselineImage, ImageType, Filter );

  tubeWrapSetMacro( BaselineNumberOfFailedPixelsTolerance, unsigned int,
    Filter );
  tubeWrapSetMacro( BaselineIntensityTolerance, PixelType, Filter );
  tubeWrapSetMacro( BaselineRadiusTolerance, unsigned int, Filter );

  tubeWrapCallMacro( ComputeBaselineDifference, Filter );

  tubeWrapGetConstObjectMacro( BaselineDifferenceImage, ImageType, Filter );
  tubeWrapGetConstObjectMacro( BaselineResampledMovingImage, ImageType,
    Filter );
  tubeWrapGetMacro( BaselineNumberOfFailedPixels, unsigned int, Filter );
  tubeWrapGetMacro( BaselineTestPassed, bool, Filter );

  //
  // Process Control
  //
  tubeWrapSetMacro( EnableLoadedRegistration, bool, Filter );
  tubeWrapGetMacro( EnableLoadedRegistration, bool, Filter );
  tubeWrapBooleanMacro( EnableLoadedRegistration, Filter );

  tubeWrapSetMacro( EnableInitialRegistration, bool, Filter );
  tubeWrapGetMacro( EnableInitialRegistration, bool, Filter );
  tubeWrapBooleanMacro( EnableInitialRegistration, Filter );

  tubeWrapSetMacro( EnableRigidRegistration, bool, Filter );
  tubeWrapGetMacro( EnableRigidRegistration, bool, Filter );
  tubeWrapBooleanMacro( EnableRigidRegistration, Filter );

  tubeWrapSetMacro( EnableAffineRegistration, bool, Filter );
  tubeWrapGetMacro( EnableAffineRegistration, bool, Filter );
  tubeWrapBooleanMacro( EnableAffineRegistration, Filter );

  tubeWrapSetMacro( EnableBSplineRegistration, bool, Filter );
  tubeWrapGetMacro( EnableBSplineRegistration, bool, Filter );
  tubeWrapBooleanMacro( EnableBSplineRegistration, Filter );

  void SetRegistration( const std::string & reg );
  void SetInterpolation( const std::string & interp );
  void SetMetric( const std::string & metric );

  //
  // Specify the Optimizer
  //
  tubeWrapSetMacro( UseEvolutionaryOptimization, bool, Filter );
  tubeWrapGetMacro( UseEvolutionaryOptimization, bool, Filter );

  tubeWrapSetMacro( ExpectedOffsetMagnitude, double, Filter );
  tubeWrapGetMacro( ExpectedOffsetMagnitude, double, Filter );

  tubeWrapSetMacro( ExpectedRotationMagnitude, double, Filter );
  tubeWrapGetMacro( ExpectedRotationMagnitude, double, Filter );

  tubeWrapSetMacro( ExpectedScaleMagnitude, double, Filter );
  tubeWrapGetMacro( ExpectedScaleMagnitude, double, Filter );

  tubeWrapSetMacro( ExpectedSkewMagnitude, double, Filter );
  tubeWrapGetMacro( ExpectedSkewMagnitude, double, Filter );

  tubeWrapSetMacro( ExpectedDeformationMagnitude, double, Filter );
  tubeWrapGetMacro( ExpectedDeformationMagnitude, double, Filter );

  //
  // Current product of the registration pipeline
  //
  tubeWrapGetConstObjectMacro( CurrentMatrixTransform, MatrixTransformType,
    Filter );
  tubeWrapGetConstObjectMacro( CurrentBSplineTransform, BSplineTransformType,
    Filter );
  tubeWrapGetConstObjectMacro( CurrentMovingImage, ImageType, Filter );
  tubeWrapGetConstObjectMacro( LoadedTransformResampledImage, ImageType,
    Filter );
  tubeWrapGetConstObjectMacro( MatrixTransformResampledImage, ImageType,
    Filter );
  tubeWrapGetConstObjectMacro( BSplineTransformResampledImage, ImageType,
    Filter );

  //
  // Metric value
  //
  tubeWrapGetMacro( FinalMetricValue, double, Filter );

  //
  // Report progress and minimize memory
  //
  tubeWrapSetMacro( ReportProgress, bool, Filter );
  tubeWrapGetMacro( ReportProgress, bool, Filter );
  tubeWrapBooleanMacro( ReportProgress, Filter );

  tubeWrapSetMacro( MinimizeMemory, bool, Filter );
  tubeWrapGetMacro( MinimizeMemory, bool, Filter );
  tubeWrapBooleanMacro( MinimizeMemory, Filter );

  // 
  // Loaded transform parameters
  //
  void LoadTransform( const std::string & transform,
    bool invertLoadedTransform=false );
  tubeWrapCallWithConstReferenceArgMacro( SaveTransform, std::string, Filter );

  tubeWrapCallWithConstReferenceArgMacro( SaveDisplacementField, std::string,
    Filter );

  void SetLoadedMatrixTransform( const MatrixTransformType & tfm, bool
    invert=false);
  tubeWrapGetConstObjectMacro( LoadedMatrixTransform, MatrixTransformType,
    Filter );

  tubeWrapForceSetConstReferenceMacro( LoadedBSplineTransform,
    BSplineTransformType, Filter );
  tubeWrapGetConstObjectMacro( LoadedBSplineTransform, BSplineTransformType,
    Filter );

  //
  // Initial Parameters
  //
  void SetInitialMethodEnum( const std::string & initialMethod );
  const std::string GetInitialMethodEnum( void );

  tubeWrapForceSetMacro( FixedLandmarks, LandmarkVectorType, Filter );
  tubeWrapForceSetMacro( MovingLandmarks, LandmarkVectorType, Filter );

  //
  // Rigid Parameters
  //
  tubeWrapSetMacro( RigidSamplingRatio, double, Filter );
  tubeWrapGetMacro( RigidSamplingRatio, double, Filter );

  tubeWrapSetMacro( RigidTargetError, double, Filter );
  tubeWrapGetMacro( RigidTargetError, double, Filter );

  tubeWrapSetMacro( RigidMaxIterations, unsigned int, Filter );
  tubeWrapGetMacro( RigidMaxIterations, unsigned int, Filter );

  void SetRigidMetricMethodEnum( const std::string & rigidMetricMethod );
  const std::string GetRigidMetricMethodEnum( void );

  void SetRigidInterpolationMethodEnum( const std::string & rigidInterpolationMethod );
  const std::string GetRigidInterpolationMethodEnum( void );

  tubeWrapGetConstObjectMacro( RigidTransform, RigidTransformType, Filter );

  tubeWrapGetMacro( RigidMetricValue, double, Filter );

  //
  // AffineParameters
  //
  tubeWrapSetMacro( AffineSamplingRatio, double, Filter );
  tubeWrapGetMacro( AffineSamplingRatio, double, Filter );

  tubeWrapSetMacro( AffineTargetError, double, Filter );
  tubeWrapGetMacro( AffineTargetError, double, Filter );

  tubeWrapSetMacro( AffineMaxIterations, unsigned int, Filter );
  tubeWrapGetMacro( AffineMaxIterations, unsigned int, Filter );

  void SetAffineMetricMethodEnum( const std::string & rigidMetricMethod );
  const std::string GetAffineMetricMethodEnum( void );

  void SetAffineInterpolationMethodEnum( const std::string & rigidInterpolationMethod );
  const std::string GetAffineInterpolationMethodEnum( void );

  tubeWrapGetConstObjectMacro( AffineTransform, AffineTransformType, Filter );

  tubeWrapGetMacro( AffineMetricValue, double, Filter );

  //
  // BSpline Parameters
  //
  tubeWrapSetMacro( BSplineSamplingRatio, double, Filter );
  tubeWrapGetMacro( BSplineSamplingRatio, double, Filter );

  tubeWrapSetMacro( BSplineTargetError, double, Filter );
  tubeWrapGetMacro( BSplineTargetError, double, Filter );

  tubeWrapSetMacro( BSplineMaxIterations, unsigned int, Filter );
  tubeWrapGetMacro( BSplineMaxIterations, unsigned int, Filter );

  tubeWrapSetMacro( BSplineControlPointPixelSpacing, double, Filter );
  tubeWrapGetMacro( BSplineControlPointPixelSpacing, double, Filter );

  void SetBSplineMetricMethodEnum( const std::string & bSplineMetricMethod);
  const std::string GetBSplineMetricMethodEnum(void);

  void SetBSplineInterpolationMethodEnum( const std::string & bSplineMetricMethod);
  const std::string GetBSplineInterpolationMethodEnum(void);

  tubeWrapGetConstObjectMacro( BSplineTransform, BSplineTransformType, Filter );
  tubeWrapGetMacro( BSplineMetricValue, double, Filter );

protected:
  RegisterImages( void );
  ~RegisterImages() {}

  void PrintSelf( std::ostream & os, itk::Indent indent ) const override;

private:
  /** itktubeRegisterImagesFilter parameters **/
  RegisterImages( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) override
    {};

  typename FilterType::Pointer  m_Filter;
};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeRegisterImages.hxx"
#endif

#endif // End !defined( __tubeRegisterImages_h )
