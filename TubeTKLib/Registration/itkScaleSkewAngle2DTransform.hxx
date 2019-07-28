/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
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
#ifndef itkScaleSkewAngle2DTransform_hxx
#define itkScaleSkewAngle2DTransform_hxx

#include "itkScaleSkewAngle2DTransform.h"
#include "itkMath.h"

namespace itk
{
// Constructor with default arguments
template<typename TParametersValueType>
ScaleSkewAngle2DTransform<TParametersValueType>
::ScaleSkewAngle2DTransform() :
  Superclass(ParametersDimension)
{
  m_UseSingleScale = false;
  m_Scale.Fill(NumericTraits<TParametersValueType>::OneValue());
  m_Skew.Fill(NumericTraits<TParametersValueType>::ZeroValue());
}

// Constructor with arguments
template<typename TParametersValueType>
ScaleSkewAngle2DTransform<TParametersValueType>
::ScaleSkewAngle2DTransform(unsigned int parametersDimension) :
  Superclass(parametersDimension)
{
  m_Scale.Fill(1.0);
  m_Skew.Fill(0.0);
}

// Constructor with arguments
template<typename TParametersValueType>
ScaleSkewAngle2DTransform<TParametersValueType>
::ScaleSkewAngle2DTransform(const MatrixType & matrix,
  const OutputVectorType & offset) :
  Superclass(matrix, offset)
{
  this->ComputeMatrixParameters();
}

// Directly set the matrix
template<typename TParametersValueType>
void
ScaleSkewAngle2DTransform<TParametersValueType>
::SetMatrix(const MatrixType & matrix)
{

  this->SetVarMatrix( matrix );
  this->ComputeOffset();
  this->ComputeMatrixParameters();
}

// Set Parameters
template<typename TParametersValueType>
void
ScaleSkewAngle2DTransform<TParametersValueType>
::SetParameters(const ParametersType & parameters)
{
  itkDebugMacro(<< "Setting parameters " << parameters);

  // Save parameters. Needed for proper operation of TransformUpdateParameters.
  if( &parameters != &(this->m_Parameters) )
    {
    this->m_Parameters = parameters;
    }
  this->SetVarAngle(parameters[0]);

  // Transfer the translation part
  TranslationType newTranslation;
  newTranslation[0] = parameters[1];
  newTranslation[1] = parameters[2];
  this->SetVarTranslation(newTranslation);

  // Matrix must be defined before translation so that offset can be computed
  // from translation
  if( m_UseSingleScale )
    {
    m_Scale[0] = parameters[3];
    m_Scale[1] = parameters[3];
    }
  else 
    {
    m_Scale[0] = parameters[3];
    m_Scale[1] = parameters[4];
    } 

  m_Skew[0] = parameters[5];
  m_Skew[1] = parameters[6];

  this->ComputeMatrix();
  this->ComputeOffset();

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();

  itkDebugMacro(<< "After setting parameters ");
}

//
// Get Parameters
//
// Parameters are ordered as:
//
// p[0] = angle
// p[1:2] = translation components
// p[3:4] = Scale
// p[5:6] = Skew {xy, yx}
//

template<typename TParametersValueType>
const typename ScaleSkewAngle2DTransform<TParametersValueType>::ParametersType
& ScaleSkewAngle2DTransform<TParametersValueType>
::GetParameters(void) const
  {
  itkDebugMacro(<< "Getting parameters ");

  this->m_Parameters[0] = this->GetAngle();

  this->m_Parameters[1] = this->GetTranslation()[0];
  this->m_Parameters[2] = this->GetTranslation()[1];

  this->m_Parameters[3] = this->GetScale()[0];
  if( m_UseSingleScale )
    {
    this->m_Parameters[4] = this->GetScale()[0];
    }
  else
    {
    this->m_Parameters[4] = this->GetScale()[1];
    }

  this->m_Parameters[5] = this->GetSkew()[0];
  this->m_Parameters[6] = this->GetSkew()[1];

  itkDebugMacro(<< "After getting parameters " << this->m_Parameters);

  return this->m_Parameters;
  }

template<typename TParametersValueType>
void
ScaleSkewAngle2DTransform<TParametersValueType>
::SetIdentity()
{
  m_Scale.Fill(NumericTraits<ScaleVectorValueType>::OneValue());
  m_Skew.Fill(NumericTraits<SkewVectorValueType>::ZeroValue());
  Superclass::SetIdentity();
}

template<typename TParametersValueType>
void
ScaleSkewAngle2DTransform<TParametersValueType>
::SetScale(const ScaleVectorType & scale)
{
  m_Scale = scale;
  this->ComputeMatrix();
}

template<typename TParametersValueType>
void
ScaleSkewAngle2DTransform<TParametersValueType>
::SetSkew(const SkewVectorType & skew)
{
  m_Skew = skew;
  this->ComputeMatrix();
}

// Compute the matrix
template<typename TParametersValueType>
void
ScaleSkewAngle2DTransform<TParametersValueType>
::ComputeMatrix(void)
{
  MatrixType rotationMatrix;
  const MatrixValueType ca = std::cos(this->GetAngle());
  const MatrixValueType sa = std::sin(this->GetAngle());
  rotationMatrix[0][0] = ca;
  rotationMatrix[0][1] = -sa;
  rotationMatrix[1][0] = sa;
  rotationMatrix[1][1] = ca;

  MatrixType scaleSkewMatrix;
  scaleSkewMatrix[0][0] = m_Scale[0];
  scaleSkewMatrix[0][1] = std::tan(m_Skew[1]) * m_Scale[1];
  scaleSkewMatrix[1][0] = std::tan(m_Skew[0]) * m_Scale[0];
  if( m_UseSingleScale )
    {
    scaleSkewMatrix[1][1] = m_Scale[0];
    }
  else
    {
    scaleSkewMatrix[1][1] = m_Scale[1];
    }

  MatrixType newMatrix = rotationMatrix * scaleSkewMatrix;
  this->SetVarMatrix(newMatrix);

  // For debugging purposes...
  this->ComputeMatrixParameters();
}

template<typename TParametersValueType>
void
ScaleSkewAngle2DTransform<TParametersValueType>
::ComputeMatrixParameters(void)
{
  vnl_matrix<TParametersValueType> matrix(2, 2);
  matrix = this->GetMatrix().GetVnlMatrix();
  vnl_svd<TParametersValueType>    svd(matrix);
  vnl_matrix<TParametersValueType> ortho(2, 2);
  ortho = svd.U() * svd.V().transpose();

  double angle = std::acos(ortho[0][0]);
  if( ortho[1][0] < 0.0 )
    {
    angle = -angle;
    }
  if( std::abs(ortho[1][0] - std::sin(angle)) > 0.000001 )
    {
    itkWarningMacro( "Bad Rotation Matrix " << this->GetMatrix() );
    }

  // Set Rotation
  //std::cout << "ComputeMatrixParameters" << std::endl;
  //std::cout << "   Original angle = " << this->GetAngle()
    //<< " : New = " << angle << std::endl;
  Superclass::SetVarAngle( angle );

  MatrixType scaleSkew = matrix * vnl_inverse(ortho);
  //std::cout << "   Original scale = " << m_Scale[0] << ", " << m_Scale[1]
    //<< " : New = " << scaleSkew[0][0] << ", " << scaleSkew[1][1] << std::endl;
  m_Scale[0] = scaleSkew[0][0];
  m_Scale[1] = scaleSkew[1][1];

  // Set Skew
  //std::cout << "   Original skew = " << m_Skew[0] << ", " << m_Skew[1]
    //<< " : New = " << std::atan(scaleSkew[1][0]/scaleSkew[0][0]) << ", "
    //<< std::atan(scaleSkew[0][1]/scaleSkew[1][1]) << std::endl;
  m_Skew[0] = std::atan2(scaleSkew[1][0], scaleSkew[0][0]);
  m_Skew[1] = std::atan2(scaleSkew[0][1], scaleSkew[1][1]);
}

// Print self
template<typename TParametersValueType>
void
ScaleSkewAngle2DTransform<TParametersValueType>
::PrintSelf(std::ostream & os, Indent indent) const
{
  Rigid2DTransform<TParametersValueType>::PrintSelf(os, indent);
  if( m_UseSingleScale )
    {
    os << indent << "UseSingleScale: true" << std::endl;
    }
  else
    {
    os << indent << "UseSingleScale: false" << std::endl;
    }
  os << indent << "Scale:          " << m_Scale        << std::endl;
  os << indent << "Skew:           " << m_Skew         << std::endl;
}

template<typename TParametersValueType>
void
ScaleSkewAngle2DTransform<TParametersValueType>
::ComputeJacobianWithRespectToParameters(const InputPointType & p, JacobianType & j) const
{
  j.SetSize( OutputSpaceDimension, this->GetNumberOfLocalParameters() );
  j.Fill(0.0);

  const double ca = std::cos( this->GetAngle() );
  const double sa = std::sin( this->GetAngle() );

  const double cx = this->GetCenter()[0];
  const double cy = this->GetCenter()[1];

  // derivatives with respect to the angle
  j[0][0] = -sa * ( p[0] - cx ) - ca * ( p[1] - cy );
  j[1][0] =  ca * ( p[0] - cx ) - sa * ( p[1] - cy );

  // compute derivatives for the translation part
  unsigned int blockOffset = 1;
  for( unsigned int dim = 0; dim < OutputSpaceDimension; dim++ )
    {
    j[dim][blockOffset + dim] = 1.0;
    }

  j[0][3] = cx;
  j[1][4] = cy;

  j[0][5] = cy;
  j[1][6] = cx;
}

} // namespace itk

#endif
