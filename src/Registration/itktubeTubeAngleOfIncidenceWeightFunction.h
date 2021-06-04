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

#ifndef __itktubeTubeAngleOfIncidenceWeightFunction_h
#define __itktubeTubeAngleOfIncidenceWeightFunction_h

#include <itkFunctionBase.h>

namespace itk
{

namespace tube
{

namespace Function
{

/** \class TubeAngleOfIncidenceWeightFunction
 *
 * \brief Weight tube points by their angle relative to a ultrasound beam.
 *
 * The tangent direction of a tube point is compared to the ultrasound beam
 * direction.  The output weight is given by:
 *
 * \f[
 * W( \mathbf{x} ) =
 *   1 - \alpha ( 1 - \cos^n \theta )
 * \f]
 *
 * Where:
 *
 * \f{eqnarray*}
 *   W( \mathbf{x} ) &=& \mbox{output weight}
 *   n             &=& \mbox{specifies angle dependence}
 *   \alpha        &=& \mbox{fractional importance of the angle of incidence}
 *   \cos \theta   &=& \mbox{angle of incidence with the tube}
 * \f}
 *
 * \warning Make sure to set the UltrasoundProbeOrigin.
 *
 */
template< class TTubePoint, class TWeight = double >
class TubeAngleOfIncidenceWeightFunction:
  public FunctionBase< TTubePoint, TWeight >
{
public:
  /** Standard class typedefs. */
  typedef TubeAngleOfIncidenceWeightFunction< TTubePoint, TWeight >
                                               Self;
  typedef FunctionBase< TTubePoint, TWeight >  Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( TubeAngleOfIncidenceWeightFunction, FunctionBase );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  typedef TWeight                           WeightType;
  typedef TTubePoint                        TubePointType;
  typedef typename TubePointType::PointType PointType;

  /** Set/Get the fractional importance of the angle of incicent an the tube
   * point weight. */
  itkSetClampMacro( FractionalImportance, double, 0, 1.0 );
  itkGetConstMacro( FractionalImportance, double );

  /** Set/Get the angle dependence. */
  itkSetMacro( AngleDependence, double );
  itkGetConstMacro( AngleDependence, double );

  /** Set/Get the ultrasound probe origin ( assuming a phase array or
   * curvilinear array transducer geometry. */
  itkSetMacro( UltrasoundProbeOrigin, PointType );
  itkGetConstReferenceMacro( UltrasoundProbeOrigin, PointType );

  WeightType Evaluate( const TubePointType & tubePoint ) const override
    {
    const PointType & position = tubePoint.GetPositionInObjectSpace();
    typename TubePointType::VectorType beam =
      ( position - m_UltrasoundProbeOrigin );
    beam.Normalize();
    typename TubePointType::VectorType tangent =
      tubePoint.GetTangentInObjectSpace();
    tangent.Normalize();
    const double dotProduct = beam * tangent;
    double cos_term =
      std::fabs( std::sqrt( 1.0 - dotProduct * dotProduct ) );
    if( m_AngleDependence != 1.0 )
      {
      cos_term = std::pow( cos_term, m_AngleDependence );
      }
    return static_cast< WeightType >( 1.0 -
      m_FractionalImportance * ( 1.0 - cos_term ) );
    }

protected:
  TubeAngleOfIncidenceWeightFunction( void ):
      m_FractionalImportance( 0.5 ),
      m_AngleDependence( 1.0 )
    {
    m_UltrasoundProbeOrigin.Fill( 0.0 );
    }
  ~TubeAngleOfIncidenceWeightFunction( void )
    {}

private:
  double    m_FractionalImportance;
  double    m_AngleDependence;
  PointType m_UltrasoundProbeOrigin;

  // purposely not implemented
  TubeAngleOfIncidenceWeightFunction( const Self & );
  // purposely not implemented
  void operator=( const Self & );
}; // End class TubeAngleOfIncidenceWeightFunction

} // End namespace Function

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeTubeAngleOfIncidenceWeightFunction_h )
