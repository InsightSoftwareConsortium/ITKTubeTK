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

#ifndef itkScaleSkewAngle2DTransform_h
#define itkScaleSkewAngle2DTransform_h

#include <iostream>
#include "itkRigid2DTransform.h"

namespace itk
{

/** \class ScaleSkewAngle2DTransform
 * \brief ScaleSkewAngle2DTransform of a vector space (e.g. space coordinates)
 *
 * This transform applies a rotation and translation & scale/skew
 * to the space
 *
 * The parameters for this transform can be set either using individual Set
 * methods or in serialized form using SetParameters() and SetFixedParameters().
 *
 * The serialization of the optimizable parameters is an array of 7 elements.
 * The first element is the angle representation of 2D rotation. The next 2
 * parameters defines the translation in each
 * dimension. The next 2 parameters defines scaling in each dimension.
 * The last 2 parameters defines the skew.
 *
 * The serialization of the fixed parameters is an array of 2 elements defining
 * the center of rotation.
 *
 * \ingroup ITKTransform
 */
template<typename TParametersValueType=double>
class ITK_TEMPLATE_EXPORT ScaleSkewAngle2DTransform :
  public Rigid2DTransform<TParametersValueType>
{
public:
  /** Standard class type alias. */
  using Self = ScaleSkewAngle2DTransform;
  using Superclass = Rigid2DTransform<TParametersValueType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ScaleSkewAngle2DTransform, Rigid2DTransform);

  /** Dimension of parameters. */
  itkStaticConstMacro(InputSpaceDimension, unsigned int, 2);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, 2);
  itkStaticConstMacro(ParametersDimension, unsigned int, 7);

  /** Parameters Type   */
  using ParametersType = typename Superclass::ParametersType;
  using FixedParametersType = typename Superclass::FixedParametersType;
  using JacobianType = typename Superclass::JacobianType;
  using ScalarType = typename Superclass::ScalarType;
  using InputPointType = typename Superclass::InputPointType;
  using OutputPointType = typename Superclass::OutputPointType;
  using InputVectorType = typename Superclass::InputVectorType;
  using OutputVectorType = typename Superclass::OutputVectorType;
  using InputVnlVectorType = typename Superclass::InputVnlVectorType;
  using OutputVnlVectorType = typename Superclass::OutputVnlVectorType;
  using InputCovariantVectorType = typename Superclass::InputCovariantVectorType;
  using OutputCovariantVectorType = typename Superclass::OutputCovariantVectorType;
  using MatrixType = typename Superclass::MatrixType;
  using MatrixValueType = typename Superclass::MatrixValueType;
  using InverseMatrixType = typename Superclass::InverseMatrixType;
  using CenterType = typename Superclass::CenterType;
  using OffsetType = typename Superclass::OffsetType;
  using TranslationType = typename Superclass::TranslationType;

  /** Scale & Skew Vector Type. */
  using ScaleVectorType = Vector<TParametersValueType, 2>;
  using SkewVectorType = Vector<TParametersValueType, 2>;

  using ScaleVectorValueType = typename ScaleVectorType::ValueType;
  using SkewVectorValueType = typename SkewVectorType::ValueType;
  using TranslationValueType = typename TranslationType::ValueType;

  using ParametersValueType = typename Superclass::ParametersValueType;

  /** Directly set the matrix of the transform.
   *
   * Orthogonality testing is bypassed in this case.
   *
   * \sa MatrixOffsetTransformBase::SetMatrix() */
  virtual void SetMatrix(const MatrixType & matrix) override;
  virtual void SetMatrix(const MatrixType & matrix, 
    const TParametersValueType tolerance) override;

  /** Set the transformation from a container of parameters
   * This is typically used by optimizers.
   * There are 15 parameters:
   *   0     angle
   *   1-2   translation
   *   3-4   Scale
   *   5-6   Skew
   **  */
  virtual void SetParameters(const ParametersType & parameters) override;

  virtual const ParametersType & GetParameters(void) const override;

  itkGetMacro( UseSingleScale, bool );
  itkSetMacro( UseSingleScale, bool );

  void SetScale(const ScaleVectorType & scale);

  itkGetConstReferenceMacro(Scale, ScaleVectorType);

  void SetSkew(const SkewVectorType & skew);

  itkGetConstReferenceMacro(Skew, SkewVectorType);

  void SetIdentity() override;

  /** This method computes the Jacobian matrix of the transformation.
   * given point or vector, returning the transformed point or
   * vector. The rank of the Jacobian will also indicate if the
   * transform is invertible at this point. */
  virtual void ComputeJacobianWithRespectToParameters( const InputPointType  & p, JacobianType & jacobian) const override;

protected:
  ScaleSkewAngle2DTransform();
  ScaleSkewAngle2DTransform(const MatrixType & matrix,
    const OutputVectorType & offset);
  ScaleSkewAngle2DTransform(unsigned int paramDims);

  ~ScaleSkewAngle2DTransform()
  {
  }

  virtual void PrintSelf(std::ostream & os, Indent indent) const override;

  void SetVarScale(const ScaleVectorType & scale)
  {
    m_Scale = scale;
  }

  void SetVarSkew(const SkewVectorType & skew)
  {
    m_Skew = skew;
  }

  /** Compute the components of the rotation matrix in the superclass. */
  void ComputeMatrix(void) override;

  void ComputeMatrixParameters(void) override;

private:
  ITK_DISALLOW_COPY_AND_MOVE(ScaleSkewAngle2DTransform);

  /**  If true, parameters[3] is used for scaling in x and y. */
  bool m_UseSingleScale;

  /**  Vector containing the scale. */
  ScaleVectorType m_Scale;

  /**  Vector containing the skew */
  SkewVectorType m_Skew;
}; // class ScaleSkewAngle2DTransform
}  // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkScaleSkewAngle2DTransform.hxx"
#endif

#endif /* __ScaleSkewAngle2DTransform_h */
