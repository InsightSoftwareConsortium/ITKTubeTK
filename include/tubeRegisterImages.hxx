/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/
#ifndef __tubeRegisterImages_hxx
#define __tubeRegisterImages_hxx

#include "tubeRegisterImages.h"

namespace tube
{

template< class TPixel, unsigned int VDimension >
RegisterImages< TPixel, VDimension >
::RegisterImages( void )
{
  m_Filter = FilterType::New();
}

template< class TPixel, unsigned int VDimension >
const typename RegisterImages<TPixel, VDimension>::ImageType *
RegisterImages< TPixel, VDimension >
::ResampleImage( std::string inter,
  const ImageType * movingImage,
  const MatrixTransformType * matrixTransform,
  //const BSplineTransformType * bsplineTransform,
  PixelType defaultPixelValue,
  double portionOfTransformToApply )
{
  const BSplineTransformType * bsplineTransform = nullptr;

  if( !strcmp(inter.c_str(), "LINEAR_INTERPOLATION" ) )
    {
    return m_Filter->ResampleImage(
      FilterType::OptimizedRegistrationMethodType::LINEAR_INTERPOLATION,
      movingImage, matrixTransform, bsplineTransform, defaultPixelValue,
      portionOfTransformToApply );
    }
  else if( !strcmp(inter.c_str(), "BSPLINE_INTERPOLATION" ) )
    {
    return m_Filter->ResampleImage(
      FilterType::OptimizedRegistrationMethodType::BSPLINE_INTERPOLATION,
      movingImage, matrixTransform, bsplineTransform, defaultPixelValue,
      portionOfTransformToApply );
    }
  else if( !strcmp(inter.c_str(), "SINC_INTERPOLATION" ) )
    {
    return m_Filter->ResampleImage(
      FilterType::OptimizedRegistrationMethodType::SINC_INTERPOLATION,
      movingImage, matrixTransform, bsplineTransform, defaultPixelValue,
      portionOfTransformToApply );
    }
  else // Assume Nearest Neighbor
    {
    return m_Filter->ResampleImage(
      FilterType::OptimizedRegistrationMethodType::NEAREST_NEIGHBOR_INTERPOLATION,
      movingImage, matrixTransform, bsplineTransform, defaultPixelValue,
      portionOfTransformToApply );
    }
}

template< class TPixel, unsigned int VDimension >
void
RegisterImages< TPixel, VDimension >
::SetRegistration( const std::string & reg )
{
  if( !strcmp(reg.c_str(), "NONE") )
    {
    m_Filter->SetRegistration( FilterType::RegistrationMethodEnumType::NONE );
    }
  else if( !strcmp(reg.c_str(), "INITIAL") )
    {
    m_Filter->SetRegistration(
      FilterType::RegistrationMethodEnumType::INITIAL );
    }
  else if( !strcmp(reg.c_str(), "RIGID") )
    {
    m_Filter->SetRegistration(
      FilterType::RegistrationMethodEnumType::RIGID );
    }
  else if( !strcmp(reg.c_str(), "AFFINE") )
    {
    m_Filter->SetRegistration(
      FilterType::RegistrationMethodEnumType::AFFINE );
    }
  else if( !strcmp(reg.c_str(), "BSPLINE") )
    {
    m_Filter->SetRegistration(
      FilterType::RegistrationMethodEnumType::BSPLINE );
    }
  else if( !strcmp(reg.c_str(), "PIPELINE_AFFINE") )
    {
    m_Filter->SetRegistration(
      FilterType::RegistrationMethodEnumType::PIPELINE_AFFINE );
    }
  else if( !strcmp(reg.c_str(), "PIPELINE_BSPLINE") )
    {
    m_Filter->SetRegistration(
      FilterType::RegistrationMethodEnumType::PIPELINE_BSPLINE );
    }
  else // if( !strcmp(reg.c_str(), "PIPELINE_RIGID") )
    {
    m_Filter->SetRegistration(
      FilterType::RegistrationMethodEnumType::PIPELINE_RIGID );
    }
}

template< class TPixel, unsigned int VDimension >
void
RegisterImages< TPixel, VDimension >
::SetInterpolation( const std::string & interp )
{
  if( !strcmp(interp.c_str(), "LINEAR_INTERPOLATIONS") )
    {
    m_Filter->SetInterpolation(
      FilterType::InterpolationMethodEnumType::LINEAR_INTERPOLATION );
    }
  else if( !strcmp(interp.c_str(), "BSPLINE_INTERPOLATION") )
    {
    m_Filter->SetInterpolation(
      FilterType::InterpolationMethodEnumType::BSPLINE_INTERPOLATION );
    }
  else if( !strcmp(interp.c_str(), "SINC_INTERPOLATION") )
    {
    m_Filter->SetInterpolation(
      FilterType::InterpolationMethodEnumType::SINC_INTERPOLATION );
    }
  else // if( !strcmp(interp.c_str(), "NEAREST_NEIGHBOR_INTERPOLATION") )
    {
    m_Filter->SetInterpolation(
      FilterType::InterpolationMethodEnumType::NEAREST_NEIGHBOR_INTERPOLATION );
    }
}

template< class TPixel, unsigned int VDimension >
void
RegisterImages< TPixel, VDimension >
::SetMetric( const std::string & metric )
{
  if( !strcmp(metric.c_str(), "MATTES_MI_METRIC") )
    {
    m_Filter->SetMetric(
      FilterType::MetricMethodEnumType::MATTES_MI_METRIC );
    }
  else if( !strcmp(metric.c_str(), "NORMALIZED_CORRELATION_METRIC") )
    {
    m_Filter->SetMetric(
      FilterType::MetricMethodEnumType::NORMALIZED_CORRELATION_METRIC );
    }
  else // if( !strcmp(metric.c_str(), "MEAN_SQUARED_ERROR_METRIC") )
    {
    m_Filter->SetMetric(
      FilterType::MetricMethodEnumType::MEAN_SQUARED_ERROR_METRIC );
    }
}

template< class TPixel, unsigned int VDimension >
const typename RegisterImages<TPixel, VDimension>::ImageType *
RegisterImages< TPixel, VDimension >
::GetFinalMovingImage( std::string inter )
{
  if( !strcmp(inter.c_str(), "LINEAR_INTERPOLATION" ) )
    {
    return m_Filter->GetFinalMovingImage( FilterType::
      OptimizedRegistrationMethodType::LINEAR_INTERPOLATION );
    }
  else if( !strcmp(inter.c_str(), "BSPLINE_INTERPOLATION" ) )
    {
    return m_Filter->GetFinalMovingImage( FilterType::
      OptimizedRegistrationMethodType::BSPLINE_INTERPOLATION );
    }
  else if( !strcmp(inter.c_str(), "SINC_INTERPOLATION" ) )
    {
    return m_Filter->GetFinalMovingImage( FilterType::
      OptimizedRegistrationMethodType::SINC_INTERPOLATION );
    }
  else // Assume Nearest Neighbor
    {
    return m_Filter->GetFinalMovingImage( FilterType::
      OptimizedRegistrationMethodType::NEAREST_NEIGHBOR_INTERPOLATION );
    }
}


template< class TPixel, unsigned int VDimension >
void
RegisterImages< TPixel, VDimension >
::LoadTransform( const std::string & t, bool invertLoadedTransform )
{
  m_Filter->LoadTransform( t, invertLoadedTransform );
}

template< class TPixel, unsigned int VDimension >
void
RegisterImages< TPixel, VDimension >
::SetLoadedMatrixTransform( const MatrixTransformType & tfm, bool invert )
{
  m_Filter->SetLoadedMatrixTransform( tfm, invert );
}

template< class TPixel, unsigned int VDimension >
void
RegisterImages< TPixel, VDimension >
::SetInitialMethodEnum( const std::string & initialMethod )
{
  if( !strcmp(initialMethod.c_str(), "INIT_WITH_CURRENT_RESULTS" ) )
  {
    m_Filter->SetInitialMethodEnum( FilterType::INIT_WITH_CURRENT_RESULTS );
  }
  else if( !strcmp(initialMethod.c_str(), "INIT_WITH_IMAGE_CENTERS" ) )
  {
    m_Filter->SetInitialMethodEnum( FilterType::INIT_WITH_IMAGE_CENTERS );
  }
  else if( !strcmp(initialMethod.c_str(), "INIT_WITH_CENTERS_OF_MASS" ) )
  {
    m_Filter->SetInitialMethodEnum( FilterType::INIT_WITH_CENTERS_OF_MASS );
  }
  else if( !strcmp(initialMethod.c_str(), "INIT_WITH_SECOND_MOMENTS" ) )
  {
    m_Filter->SetInitialMethodEnum( FilterType::INIT_WITH_SECOND_MOMENTS );
  }
  else if( !strcmp(initialMethod.c_str(), "INIT_WITH_LANDMARKS" ) )
  {
    m_Filter->SetInitialMethodEnum( FilterType::INIT_WITH_LANDMARKS );
  }
  else // Default: if( !strcmp(inter.c_str(), "INIT_WITH_NONE" ) )
  {
    m_Filter->SetInitialMethodEnum( FilterType::INIT_WITH_NONE );
  }
}

template< class TPixel, unsigned int VDimension >
const std::string
RegisterImages< TPixel, VDimension >
::GetInitialMethodEnum( void )
{
  if( m_Filter->GetInitialMethodEnum() ==
      FilterType::INIT_WITH_CURRENT_RESULTS )
  {
    return "INIT_WITH_CURRENT_RESULTS";
  }
  else if( m_Filter->GetInitialMethodEnum() ==
           FilterType::INIT_WITH_IMAGE_CENTERS )
  {
    return "INIT_WITH_IMAGE_CENTERS";
  }
  else if( m_Filter->GetInitialMethodEnum() ==
           FilterType::INIT_WITH_CENTERS_OF_MASS )
  {
    return "INIT_WITH_CENTERS_OF_MASS";
  }
  else if( m_Filter->GetInitialMethodEnum() ==
           FilterType::INIT_WITH_SECOND_MOMENTS )
  {
    return "INIT_WITH_SECOND_MOMENTS";
  }
  else if( m_Filter->GetInitialMethodEnum() ==
           FilterType::INIT_WITH_LANDMARKS )
  {
    return "INIT_WITH_LANDMARKS";
  }
  else // Default: if( !strcmp(inter.c_str(), "INIT_WITH_NONE" ) )
  {
    return "INIT_WITH_NONE";
  }
}

template< class TPixel, unsigned int VDimension >
void
RegisterImages< TPixel, VDimension >
::SetRigidMetricMethodEnum( const std::string & rigidMetricMethod )
{
  if( !strcmp(rigidMetricMethod.c_str(), "NORMALIZED_CORRELATION_METRIC" ) )
  {
    m_Filter->SetRigidMetricMethodEnum( 
      FilterType::MetricMethodEnumType::NORMALIZED_CORRELATION_METRIC );
  }
  else if( !strcmp(rigidMetricMethod.c_str(), "MEAN_SQUARED_ERROR_METRIC" ) )
  {
    m_Filter->SetRigidMetricMethodEnum( 
      FilterType::MetricMethodEnumType::MEAN_SQUARED_ERROR_METRIC );
  }
  else // Default: if( !strcmp(rigidMetricMethod.c_str(), "MATTES_MI_METRIC" ) )
  {
    m_Filter->SetRigidMetricMethodEnum( 
      FilterType::MetricMethodEnumType::MATTES_MI_METRIC );
  }
}

template< class TPixel, unsigned int VDimension >
const std::string
RegisterImages< TPixel, VDimension >
::GetRigidMetricMethodEnum( void )
{
  if( m_Filter->GetRigidMetricMethodEnum() == 
    FilterType::MetricMethodEnumType::NORMALIZED_CORRELATION_METRIC )
  {
    return "NORMALIZED_CORRELATION_METRIC";
  }
  else if( m_Filter->GetRigidMetricMethodEnum() == 
    FilterType::MetricMethodEnumType::MEAN_SQUARED_ERROR_METRIC )
  {
    return "MEAN_SQUARED_ERROR_METRIC";
  }
  else // Default: if( !strcmp(inter.c_str(), "MATTES_MI_METRIC" ) )
  {
    return "MATTES_MI_METRIC";
  }
}

template< class TPixel, unsigned int VDimension >
void
RegisterImages< TPixel, VDimension >
::SetRigidInterpolationMethodEnum( const std::string & rigidInterpolationMethod )
{
  if( !strcmp(rigidInterpolationMethod.c_str(), "LINEAR_INTERPOLATION" ) )
  {
    m_Filter->SetRigidInterpolationMethodEnum( 
      FilterType::InterpolationMethodEnumType::LINEAR_INTERPOLATION );
  }
  else if( !strcmp(rigidInterpolationMethod.c_str(), "BSPLINE_INTERPOLATION" ) )
  {
    m_Filter->SetRigidInterpolationMethodEnum( 
      FilterType::InterpolationMethodEnumType::BSPLINE_INTERPOLATION );
  }
  else if( !strcmp(rigidInterpolationMethod.c_str(), "SINC_INTERPOLATION" ) )
  {
    m_Filter->SetRigidInterpolationMethodEnum( 
      FilterType::InterpolationMethodEnumType::SINC_INTERPOLATION );
  }
  else // Default: if( !strcmp(rigidInterpolationMethod.c_str(),
       //   "NEAREST_NEIGHBOR_INTERPOLATION" ) )
  {
    m_Filter->SetRigidInterpolationMethodEnum( 
      FilterType::InterpolationMethodEnumType::NEAREST_NEIGHBOR_INTERPOLATION );
  }
}

template< class TPixel, unsigned int VDimension >
const std::string
RegisterImages< TPixel, VDimension >
::GetRigidInterpolationMethodEnum( void )
{
  if( m_Filter->GetRigidInterpolationMethodEnum() == 
    FilterType::InterpolationMethodEnumType::LINEAR_INTERPOLATION )
  {
    return "LINEAR_INTERPOLATION";
  }
  else if( m_Filter->GetRigidInterpolationMethodEnum() == 
    FilterType::InterpolationMethodEnumType::BSPLINE_INTERPOLATION )
  {
    return "BSPLINE_INTERPOLATION";
  }
  else if( m_Filter->GetRigidInterpolationMethodEnum() == 
    FilterType::InterpolationMethodEnumType::SINC_INTERPOLATION )
  {
    return "SINC_INTERPOLATION";
  }
  else // Default
  {
    return "NEAREST_NEIGHBOR_INTERPOLATION";
  }
}

template< class TPixel, unsigned int VDimension >
void
RegisterImages< TPixel, VDimension >
::SetAffineMetricMethodEnum( const std::string & affineMetricMethod )
{
  if( !strcmp(affineMetricMethod.c_str(), "NORMALIZED_CORRELATION_METRIC" ) )
  {
    m_Filter->SetAffineMetricMethodEnum( 
      FilterType::MetricMethodEnumType::NORMALIZED_CORRELATION_METRIC );
  }
  else if( !strcmp(affineMetricMethod.c_str(), "MEAN_SQUARED_ERROR_METRIC" ) )
  {
    m_Filter->SetAffineMetricMethodEnum( 
      FilterType::MetricMethodEnumType::MEAN_SQUARED_ERROR_METRIC );
  }
  else // Default: if( !strcmp(affineMetricMethod.c_str(),
       //   "MATTES_MI_METRIC" ) )
  {
    m_Filter->SetAffineMetricMethodEnum( 
      FilterType::MetricMethodEnumType::MATTES_MI_METRIC );
  }
}

template< class TPixel, unsigned int VDimension >
const std::string
RegisterImages< TPixel, VDimension >
::GetAffineMetricMethodEnum( void )
{
  if( m_Filter->GetAffineMetricMethodEnum() == 
    FilterType::MetricMethodEnumType::NORMALIZED_CORRELATION_METRIC )
  {
    return "NORMALIZED_CORRELATION_METRIC";
  }
  else if( m_Filter->GetAffineMetricMethodEnum() == 
    FilterType::MetricMethodEnumType::MEAN_SQUARED_ERROR_METRIC )
  {
    return "MEAN_SQUARED_ERROR_METRIC";
  }
  else // Default: if( !strcmp(inter.c_str(), "MATTES_MI_METRIC" ) )
  {
    return "MATTES_MI_METRIC";
  }
}

template< class TPixel, unsigned int VDimension >
void 
RegisterImages< TPixel, VDimension >
::SetAffineInterpolationMethodEnum(
  const std::string & affineInterpolationMethod )
{
  if( !strcmp(affineInterpolationMethod.c_str(), "LINEAR_INTERPOLATION" ) )
  {
    m_Filter->SetAffineInterpolationMethodEnum( 
      FilterType::InterpolationMethodEnumType::LINEAR_INTERPOLATION );
  }
  else if( !strcmp(affineInterpolationMethod.c_str(), "BSPLINE_INTERPOLATION" ) )
  {
    m_Filter->SetAffineInterpolationMethodEnum( 
      FilterType::InterpolationMethodEnumType::BSPLINE_INTERPOLATION );
  }
  else if( !strcmp(affineInterpolationMethod.c_str(), "SINC_INTERPOLATION" ) )
  {
    m_Filter->SetAffineInterpolationMethodEnum( 
      FilterType::InterpolationMethodEnumType::SINC_INTERPOLATION );
  }
  else // Default: if( !strcmp(affineInterpolationMethod.c_str(),
       //   "NEAREST_NEIGHBOR_INTERPOLATION" ) )
  {
    m_Filter->SetAffineInterpolationMethodEnum( 
      FilterType::InterpolationMethodEnumType::NEAREST_NEIGHBOR_INTERPOLATION );
  }
}

template< class TPixel, unsigned int VDimension >
const std::string
RegisterImages< TPixel, VDimension >
::GetAffineInterpolationMethodEnum( void )
{
  if( m_Filter->GetAffineInterpolationMethodEnum() == 
    FilterType::InterpolationMethodEnumType::LINEAR_INTERPOLATION )
  {
    return "LINEAR_INTERPOLATION";
  }
  else if( m_Filter->GetAffineInterpolationMethodEnum() == 
    FilterType::InterpolationMethodEnumType::BSPLINE_INTERPOLATION )
  {
    return "BSPLINE_INTERPOLATION";
  }
  else if( m_Filter->GetAffineInterpolationMethodEnum() == 
    FilterType::InterpolationMethodEnumType::SINC_INTERPOLATION )
  {
    return "SINC_INTERPOLATION";
  }
  else // Default
  {
    return "NEAREST_NEIGHBOR_INTERPOLATION";
  }
}

template< class TPixel, unsigned int VDimension >
void
RegisterImages< TPixel, VDimension >
::SetBSplineMetricMethodEnum( const std::string & bSplineMetricMethod )
{
  if( !strcmp(bSplineMetricMethod.c_str(), "NORMALIZED_CORRELATION_METRIC" ) )
  {
    m_Filter->SetBSplineMetricMethodEnum( 
      FilterType::MetricMethodEnumType::NORMALIZED_CORRELATION_METRIC );
  }
  else if( !strcmp(bSplineMetricMethod.c_str(), "MEAN_SQUARED_ERROR_METRIC" ) )
  {
    m_Filter->SetBSplineMetricMethodEnum( 
      FilterType::MetricMethodEnumType::MEAN_SQUARED_ERROR_METRIC );
  }
  else // Default: if( !strcmp(bSplineMetricMethod.c_str(),
       //   "MATTES_MI_METRIC" ) )
  {
    m_Filter->SetBSplineMetricMethodEnum( 
      FilterType::MetricMethodEnumType::MATTES_MI_METRIC );
  }
}

template< class TPixel, unsigned int VDimension >
const std::string
RegisterImages< TPixel, VDimension >
::GetBSplineMetricMethodEnum( void )
{
  if( m_Filter->GetBSplineMetricMethodEnum() == 
    FilterType::MetricMethodEnumType::NORMALIZED_CORRELATION_METRIC )
  {
    return "NORMALIZED_CORRELATION_METRIC";
  }
  else if( m_Filter->GetBSplineMetricMethodEnum() == 
    FilterType::MetricMethodEnumType::MEAN_SQUARED_ERROR_METRIC )
  {
    return "MEAN_SQUARED_ERROR_METRIC";
  }
  else // Default: if( !strcmp(inter.c_str(), "MATTES_MI_METRIC" ) )
  {
    return "MATTES_MI_METRIC";
  }
}

template< class TPixel, unsigned int VDimension >
void 
RegisterImages< TPixel, VDimension >
::SetBSplineInterpolationMethodEnum(
  const std::string & bSplineInterpolationMethod )
{
  if( !strcmp(bSplineInterpolationMethod.c_str(), "LINEAR_INTERPOLATION" ) )
  {
    m_Filter->SetBSplineInterpolationMethodEnum( 
      FilterType::InterpolationMethodEnumType::LINEAR_INTERPOLATION );
  }
  else if( !strcmp(bSplineInterpolationMethod.c_str(), "BSPLINE_INTERPOLATION" ) )
  {
    m_Filter->SetBSplineInterpolationMethodEnum( 
      FilterType::InterpolationMethodEnumType::BSPLINE_INTERPOLATION );
  }
  else if( !strcmp(bSplineInterpolationMethod.c_str(), "SINC_INTERPOLATION" ) )
  {
    m_Filter->SetBSplineInterpolationMethodEnum( 
      FilterType::InterpolationMethodEnumType::SINC_INTERPOLATION );
  }
  else // Default: if( !strcmp(bSplineInterpolationMethod.c_str(),
       //   "NEAREST_NEIGHBOR_INTERPOLATION" ) )
  {
    m_Filter->SetBSplineInterpolationMethodEnum( 
      FilterType::InterpolationMethodEnumType::NEAREST_NEIGHBOR_INTERPOLATION );
  }
}

template< class TPixel, unsigned int VDimension >
const std::string
RegisterImages< TPixel, VDimension >
::GetBSplineInterpolationMethodEnum( void )
{
  if( m_Filter->GetBSplineInterpolationMethodEnum() == 
    FilterType::InterpolationMethodEnumType::LINEAR_INTERPOLATION )
  {
    return "LINEAR_INTERPOLATION";
  }
  else if( m_Filter->GetBSplineInterpolationMethodEnum() == 
    FilterType::InterpolationMethodEnumType::BSPLINE_INTERPOLATION )
  {
    return "BSPLINE_INTERPOLATION";
  }
  else if( m_Filter->GetBSplineInterpolationMethodEnum() == 
    FilterType::InterpolationMethodEnumType::SINC_INTERPOLATION )
  {
    return "SINC_INTERPOLATION";
  }
  else // Default
  {
    return "NEAREST_NEIGHBOR_INTERPOLATION";
  }
}

template< class TPixel, unsigned int VDimension >
void
RegisterImages< TPixel, VDimension >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << m_Filter << std::endl;
}

}; // end namespace tube

#endif
