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

#ifndef __itkUltrasoundProbeGeometryCalculator_h
#define __itkUltrasoundProbeGeometryCalculator_h

#include <itkProcessObject.h>
#include <itkSimpleDataObjectDecorator.h>

namespace itk
{

namespace tube
{

/** \class UltrasoundProbeGeometryCalculator
 *
 * \brief Compute probe origin and radius of curvature from a scan-converted
 * image.
 *
 * Given a scan-converted ultrasound B-Mode "sector" image that was
 * generated
 * with a curvilinear or phased array probe, compute a center of curvature
 * ( UltrasoundProbeOrigin ) and the distance to the start of acquisition
 * ( StartOfAcquisitionRadius ).
 *
 * We iterate in from the sides of the image orthogonal to the
 * GeneralBeamDirection and find the first values that are not equal to the
 * BackgroundValue.  The parameters of lines for both sides of the image are
 * created for each sector edge point.  The median values of these line
 * parameters define two lines in one plane.  The output
 * UltrasoundProbeOrigin
 * is the intersection of these lines.  By radiating out from this origin,
 * the
 * median of the first non-background values determine the output
 * StartOfAcquisitionRadius.
 */
template< class TInputImage >
class UltrasoundProbeGeometryCalculator : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef UltrasoundProbeGeometryCalculator  Self;
  typedef ProcessObject                      Superclass;
  typedef SmartPointer< Self >               Pointer;
  typedef SmartPointer< const Self >         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( UltrasoundProbeGeometryCalculator, ProcessObject );

  /** Some convenient typedefs. */
  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::PixelType      InputPixelType;
  typedef typename InputImageType::PointType      OriginType;
  typedef double                                  RadiusType;

  typedef SimpleDataObjectDecorator< OriginType > DecoratedOriginType;
  typedef SimpleDataObjectDecorator< RadiusType > DecoratedRadiusType;

  /** Set/Get the general beam direction.  This is the approximate direction
   * that the beam is orientated in the sector image. Defaults to zero.  */
  itkSetMacro( GeneralBeamDirection, unsigned int );
  itkGetConstMacro( GeneralBeamDirection, unsigned int );

  /** Set/Get the background value that is used in the input image outside
   * of the scanned sector. */
  itkSetMacro( BackgroundValue, InputPixelType );
  itkGetConstMacro( BackgroundValue, InputPixelType );

  virtual void SetInput( const InputImageType * image );
  virtual const InputImageType * GetInput( void ) const;

  const OriginType & GetUltrasoundProbeOrigin( void ) const;
  const RadiusType & GetStartOfAcquisitionRadius( void ) const;

protected:
  UltrasoundProbeGeometryCalculator( void );
  virtual ~UltrasoundProbeGeometryCalculator( void );

  virtual void PrintSelf( std::ostream & os, Indent indent ) const override;

  using Superclass::MakeOutput;
  virtual DataObject::Pointer MakeOutput(
    DataObjectPointerArraySizeType index ) override;

  virtual void GenerateData( void ) override;

private:
  // purposely not implemented
  UltrasoundProbeGeometryCalculator( const Self & );
  // purposely not implemented
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) override
    {};

  unsigned int   m_GeneralBeamDirection;
  InputPixelType m_BackgroundValue;

}; // End class UltrasoundProbeGeometryCalculator

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkUltrasoundProbeGeometryCalculator.hxx"
#endif

#endif // End !defined( __itkUltrasoundProbeGeometryCalculator_h )
