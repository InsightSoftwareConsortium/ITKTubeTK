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

#ifndef __itkGradientBasedAngleOfIncidenceImageFilter_h
#define __itkGradientBasedAngleOfIncidenceImageFilter_h

#include <itkCastImageFilter.h>
#include <itkCovariantVector.h>
#include <itkImageToImageFilter.h>

namespace itk
{
/** \class GradientBasedAngleOfIncidenceImageFilter
 * \brief Computes cosine of the angle of incidence.
 *
 * The cosine of the angle of incidence is computed as the angle between the
 * ultrasound beam direction and the normal of the "local surface", which is
 * computed as the local gradient.
 *
 * For every input pixel location, the beam direction is computed by normalizing
 * the vector from that location to the center of rotation of a phased array or
 * curvilinear array probe, specified with the \c UltrasoundProbeOrigin.  The
 * gradient is computed with a gradient filter of the user's choice -- the
 * default is a GradientImageFilter, but a difference filter, e.g. a
 * GradientRecursiveGaussianImageFilter could be used instead.
 *
 * The cosine of the angle of incidence is computed as the dot product of the
 * two normalized vectors.
 *
 * \ingroup ImageToImageFilter
 */
template< class TInputImage, class TOutputImage, class TOperatorValue = float >
class GradientBasedAngleOfIncidenceImageFilter
  : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef GradientBasedAngleOfIncidenceImageFilter            Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage >     Superclass;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( GradientBasedAngleOfIncidenceImageFilter, ImageToImageFilter );

  /** Some convenient typedefs. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::PointType    OriginType;
  typedef TOutputImage                          OutputImageType;
  typedef typename OutputImageType::RegionType  OutputImageRegionType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    InputImageType::ImageDimension );

  typedef TOperatorValue
    OperatorValueType;
  typedef Image< OperatorValueType, ImageDimension >
    OperatorImageType;
  typedef CovariantVector< OperatorValueType, ImageDimension >
    GradientOutputPixelType;
  typedef Image< GradientOutputPixelType, ImageDimension >
    GradientOutputImageType;
  typedef ImageToImageFilter< OperatorImageType, GradientOutputImageType >
    GradientFilterType;
  typedef Vector< TOperatorValue, InputImageType::ImageDimension >
    BeamDirectionType;

  /** Probe type.  Determines how the beam angle is calculated.  For
   * CURVILINEAR or PHASED, the UltrasoundProbeOrigin must be set.  For
   * a LINEAR probe, the UltrasoundProbeDirection must be set. */
  typedef enum
    {
    CURVILINEAR,
    PHASED,
    LINEAR
    }
  ProbeType;

  /** Set/Get the probe type.  This determines how the beam direction is
   * computed. */
  itkSetMacro( UltrasoundProbeType, ProbeType );
  itkGetConstMacro( UltrasoundProbeType, ProbeType );

  /** Set/Get the location of the ultrasound beam probe center of rotation.
   * This is only valid when the UltrasoundProbeType is CURVILINEAR or PHASED. */
  itkSetMacro( UltrasoundProbeOrigin, OriginType );
  itkGetConstMacro( UltrasoundProbeOrigin, OriginType );

  /** Set/Get the direction of the ultrasound beam.  This is only valid when the
   * UltrasoundProbeType is LINEAR. */
  void SetUltrasoundProbeBeamDirection( const BeamDirectionType & beamDirection );
  itkGetConstMacro( UltrasoundProbeBeamDirection, BeamDirectionType );

  /** Set/Get the filter used to calculate the gradients of the input image.
   * The default is a simple GradientImageFilter. */
  itkSetObjectMacro( GradientFilter, GradientFilterType );
  itkGetObjectMacro( GradientFilter, GradientFilterType );

  /** Set/Get the tolerance for the gradient magnitude.  If the gradient
   * magnitude is below this value, the output is set to zero. */
  itkSetMacro( GradientMagnitudeTolerance, double );
  itkGetConstMacro( GradientMagnitudeTolerance, double );

protected:
  GradientBasedAngleOfIncidenceImageFilter( void );
  virtual ~GradientBasedAngleOfIncidenceImageFilter( void ) {}

  virtual void PrintSelf( std::ostream & os, Indent indent ) const;

  virtual void BeforeThreadedGenerateData( void );
  virtual void ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread,
    ThreadIdType threadId );

private:
  GradientBasedAngleOfIncidenceImageFilter( const Self & ); //purposely not implemented
  void operator=( const Self & );          //purposely not implemented

  typedef CastImageFilter< InputImageType, OperatorImageType > CastImageFilterType;
  typename CastImageFilterType::Pointer m_CastImageFilter;


  typename GradientFilterType::Pointer m_GradientFilter;

  double            m_GradientMagnitudeTolerance;
  ProbeType         m_UltrasoundProbeType;
  OriginType        m_UltrasoundProbeOrigin;
  BeamDirectionType m_UltrasoundProbeBeamDirection;

}; // End class GradientBasedAngleOfIncidenceImageFilter

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGradientBasedAngleOfIncidenceImageFilter.hxx"
#endif

#endif // End !defined( __itkGradientBasedAngleOfIncidenceImageFilter_h )
