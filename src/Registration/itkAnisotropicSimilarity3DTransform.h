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

#ifndef __itkAnisotropicSimilarity3DTransform_h
#define __itkAnisotropicSimilarity3DTransform_h

#include <iostream>
#include "itkVersorRigid3DTransform.h"

namespace itk
{

/** \brief AnisotropicSimilarity3DTransform of a vector space ( e.g. space
 * coordinates )
 *
 * This transform applies a rotation, translation and anisotropic scaling
 * to the space.
 *
 * The parameters for this transform can be set either using individual
 * Set methods or in serialized form using SetParameters() and
 * SetFixedParameters().
 *
 * The serialization of the optimizable parameters is an array of 9
 * elements. The first 3 elements are the components of the versor
 * representation of 3D rotation. The next 3 parameters defines the
 * translation in each dimension. The last parameter defines the
 * anisotropic scaling.
 *
 * The serialization of the fixed parameters is an array of 3 elements
 * defining the center of rotation.
 *
 * \ingroup Transforms
 *
 * \sa VersorRigid3DTransform
 */
template <class TScalarType = double>
// Data type for scalars ( float or double )
class AnisotropicSimilarity3DTransform :
  public VersorRigid3DTransform<TScalarType>
{
public:
  /** Standard class typedefs. */
  typedef AnisotropicSimilarity3DTransform    Self;
  typedef VersorRigid3DTransform<TScalarType> Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( AnisotropicSimilarity3DTransform, VersorRigid3DTransform );

  /** Dimension of parameters. */
  itkStaticConstMacro( SpaceDimension, unsigned int, 3 );
  itkStaticConstMacro( InputSpaceDimension, unsigned int, 3 );
  itkStaticConstMacro( OutputSpaceDimension, unsigned int, 3 );
  itkStaticConstMacro( ParametersDimension, unsigned int, 9 );

  /** Parameters Type   */
  typedef typename Superclass::ParametersType      ParametersType;
  typedef typename Superclass::JacobianType        JacobianType;
  typedef typename Superclass::ScalarType          ScalarType;
  typedef typename Superclass::InputPointType      InputPointType;
  typedef typename Superclass::OutputPointType     OutputPointType;
  typedef typename Superclass::InputVectorType     InputVectorType;
  typedef typename Superclass::OutputVectorType    OutputVectorType;
  typedef typename Superclass::InputVnlVectorType  InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType OutputVnlVectorType;

  typedef typename Superclass::InputCovariantVectorType
                                                 InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType
                                                 OutputCovariantVectorType;
  typedef typename Superclass::MatrixType        MatrixType;
  typedef typename Superclass::InverseMatrixType InverseMatrixType;
  typedef typename Superclass::CenterType        CenterType;
  typedef typename Superclass::OffsetType        OffsetType;
  typedef typename Superclass::TranslationType   TranslationType;

  /** Versor type. */
  typedef typename Superclass::VersorType      VersorType;
  typedef typename Superclass::AxisType        AxisType;
  typedef typename Superclass::AngleType       AngleType;
  typedef typename Superclass::InputVectorType VectorType;
  typedef          TScalarType                 ScaleType;

  /** Directly set the rotation matrix of the transform.
   * \warning The input matrix must be orthogonal with isotropic scaling
   * to within a specified tolerance, else an exception is thrown.
   *
   * \sa MatrixOffsetTransformBase::SetMatrix() */
  virtual void SetMatrix( const MatrixType & matrix ) override;
  virtual void SetMatrix( const MatrixType & matrix, const double tolerance )
    override;

  /** Set the transformation from a container of parameters This is
   * typically used by optimizers.  There are 7 parameters. The first
   * three represent the versor, the next three represent the translation
   * and the last one represents the scaling factor. */
  void SetParameters( const ParametersType & parameters ) override;

  virtual const ParametersType & GetParameters( void ) const override;

  /** Set/Get the value of the isotropic scaling factor */
  void SetScale( ScaleType scale );

  void SetScale( VectorType scale );

  itkGetConstReferenceMacro( Scale, VectorType );

  /** This method computes the Jacobian matrix of the transformation.
   * given point or vector, returning the transformed point or
   * vector. The rank of the Jacobian will also indicate if the
   * transform is invertible at this point. */
  virtual const JacobianType & GetJacobian( const InputPointType  & point )
    const;

  virtual void ComputeJacobianWithRespectToParameters(
    const InputPointType & p, JacobianType & jacobian ) const override;

protected:
  AnisotropicSimilarity3DTransform( const MatrixType & matrix,
    const OutputVectorType & offset );
  AnisotropicSimilarity3DTransform( unsigned int paramDim );
  AnisotropicSimilarity3DTransform();
  ~AnisotropicSimilarity3DTransform()
    {
    };

  void PrintSelf( std::ostream & os, Indent indent ) const override;

  /** Recomputes the matrix by calling the Superclass::ComputeMatrix() and
   * then applying the scale factor. */
  void ComputeMatrix() override;

  /** Computes the parameters from an input matrix. */
  void ComputeMatrixParameters() override;

private:
  // purposely not implemented
  AnisotropicSimilarity3DTransform( const Self & );
  // purposely not implemented
  void operator=( const Self & );

  VectorType           m_Scale;
  mutable JacobianType m_NonThreadsafeSharedJacobian;

}; // class AnisotropicSimilarity3DTransform

}  // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAnisotropicSimilarity3DTransform.hxx"
#endif

#endif /* __itkAnisotropicSimilarity3DTransform_h */
