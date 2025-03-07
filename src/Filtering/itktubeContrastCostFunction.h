/*=========================================================================

Library:   TubeTKLib

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

#ifndef itktubeContrastCostFunction_h
#define itktubeContrastCostFunction_h

#include <itkSmoothingRecursiveGaussianImageFilter.h>
#include <itkSingleValuedCostFunction.h>
namespace itk
{

namespace tube
{

/** \class ContrastCostFunction
 */

template< class TPixel, unsigned int VDimension >
class ContrastCostFunction : public itk::SingleValuedCostFunction
{
public:

  /** Tube class type alias */
  using Self = ContrastCostFunction;
  using Superclass = SingleValuedCostFunction;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ContrastCostFunction, SingleValuedCostFunction );

  using MeasureType = Superclass::MeasureType;
  using ParametersType = Superclass::ParametersType;
  using DerivativeType = Superclass::DerivativeType;
  using ImageType = itk::Image< TPixel, VDimension >;

  using BlurFilterType = itk::SmoothingRecursiveGaussianImageFilter< ImageType, ImageType >;
  /** Set/Get input image */
  itkSetConstObjectMacro( InputImage, ImageType );
  itkGetConstObjectMacro( InputImage, ImageType );

  /** Set/Get input mask image */
  itkSetObjectMacro( InputMask, ImageType );
  itkGetModifiableObjectMacro( InputMask, ImageType );

  /** Set/Get Mask Object Value */
  itkSetMacro( MaskObjectValue, int );
  itkGetMacro( MaskObjectValue, int );

  /** Set/Get Mask Background Value */
  itkSetMacro( MaskBackgroundValue, int );
  itkGetMacro( MaskBackgroundValue, int );

  /** Set/Get output  image */
  itkSetObjectMacro( OutputImage, ImageType );
  itkGetModifiableObjectMacro( OutputImage, ImageType );

  /** Set Scales */
  void SetScales( ParametersType & scales );

  unsigned int GetNumberOfParameters( void ) const override;

  /** This method returns the value of the cost function corresponding
    * to the specified parameters.    */
  double GetValue( const ParametersType & parameters ) const override;

   /** This method returns the derivative of the cost function corresponding
    * to the specified parameters.   */
  void GetDerivative( const ParametersType & parameters,
                             DerivativeType & derivative ) const override;

  void Initialize( void );
protected:

  ContrastCostFunction( void );
  ~ContrastCostFunction( void ) {};

private:
  ContrastCostFunction( const Self & );
  void operator=( const Self & );

  typename ImageType::ConstPointer    m_InputImage;
  typename ImageType::Pointer         m_InputMask;
  mutable typename ImageType::Pointer m_OutputImage;

  double                              m_InputMean;
  int                                 m_MaskObjectValue;
  int                                 m_MaskBackgroundValue;
  ParametersType                      m_Scales;
  mutable unsigned int                m_CallsToGetValue;

}; // End class ContrastCostFunction

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeContrastCostFunction.hxx"
#endif

#endif // End !defined( __itktubeContrastCostFunction_h )
