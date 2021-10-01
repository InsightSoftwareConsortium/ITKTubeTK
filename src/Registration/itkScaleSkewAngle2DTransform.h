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
  /** Standard class typedefs. */
  typedef ScaleSkewAngle2DTransform                    Self;
  typedef Rigid2DTransform<TParametersValueType>       Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ScaleSkewAngle2DTransform, Rigid2DTransform);

  /** Dimension of parameters. */
  itkStaticConstMacro(InputSpaceDimension, unsigned int, 2);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, 2);
  itkStaticConstMacro(ParametersDimension, unsigned int, 7);

  /** Parameters Type   */
  typedef typename Superclass::ParametersType            ParametersType;
  typedef typename Superclass::FixedParametersType       FixedParametersType;
  typedef typename Superclass::JacobianType              JacobianType;
  typedef typename Superclass::ScalarType                ScalarType;
  typedef typename Superclass::InputPointType            InputPointType;
  typedef typename Superclass::OutputPointType           OutputPointType;
  typedef typename Superclass::InputVectorType           InputVectorType;
  typedef typename Superclass::OutputVectorType          OutputVectorType;
  typedef typename Superclass::InputVnlVectorType        InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType       OutputVnlVectorType;
  typedef typename Superclass::InputCovariantVectorType  InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType OutputCovariantVectorType;
  typedef typename Superclass::MatrixType                MatrixType;
  typedef typename Superclass::MatrixValueType           MatrixValueType;
  typedef typename Superclass::InverseMatrixType         InverseMatrixType;
  typedef typename Superclass::CenterType                CenterType;
  typedef typename Superclass::OffsetType                OffsetType;
  typedef typename Superclass::TranslationType           TranslationType;

  /** Scale & Skew Vector Type. */
  typedef Vector<TParametersValueType, 2> ScaleVectorType;
  typedef Vector<TParametersValueType, 2> SkewVectorType;

  typedef typename ScaleVectorType::ValueType ScaleVectorValueType;
  typedef typename SkewVectorType::ValueType  SkewVectorValueType;
  typedef typename TranslationType::ValueType TranslationValueType;

  typedef typename Superclass::ParametersValueType ParametersValueType;

  /** Directly set the matrix of the transform.
   *
   * Orthogonality testing is bypassed in this case.
   *
   * \sa MatrixOffsetTransformBase::SetMatrix() */
  virtual void SetMatrix(const MatrixType & matrix) ITK_OVERRIDE;
  virtual void SetMatrix(const MatrixType & matrix, 
    const TParametersValueType tolerance) ITK_OVERRIDE;

  /** Set the transformation from a container of parameters
   * This is typically used by optimizers.
   * There are 15 parameters:
   *   0     angle
   *   1-2   translation
   *   3-4   Scale
   *   5-6   Skew
   **  */
  virtual void SetParameters(const ParametersType & parameters) ITK_OVERRIDE;

  virtual const ParametersType & GetParameters(void) const ITK_OVERRIDE;

  itkGetMacro( UseSingleScale, bool );
  itkSetMacro( UseSingleScale, bool );

  void SetScale(const ScaleVectorType & scale);

  itkGetConstReferenceMacro(Scale, ScaleVectorType);

  void SetSkew(const SkewVectorType & skew);

  itkGetConstReferenceMacro(Skew, SkewVectorType);

  void SetIdentity() ITK_OVERRIDE;

  /** This method computes the Jacobian matrix of the transformation.
   * given point or vector, returning the transformed point or
   * vector. The rank of the Jacobian will also indicate if the
   * transform is invertible at this point. */
  virtual void ComputeJacobianWithRespectToParameters( const InputPointType  & p, JacobianType & jacobian) const ITK_OVERRIDE;

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
  void ComputeMatrix(void) ITK_OVERRIDE;

  void ComputeMatrixParameters(void) ITK_OVERRIDE;

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
