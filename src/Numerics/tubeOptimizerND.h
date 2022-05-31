/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#ifndef __tubeOptimizerND_h
#define __tubeOptimizerND_h

#include "tubeMatrixMath.h"
#include "tubeOptimizer1D.h"
#include "tubeUserFunction.h"

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

namespace tube
{

class OptimizerND : public Object
{
public:

  typedef OptimizerND                             Self;
  typedef Object                                  Superclass;
  typedef Self *                                  Pointer;
  typedef const Self *                            ConstPointer;

  typedef vnl_matrix< double >                    MatrixType;
  typedef vnl_vector< double >                    VectorType;
  typedef UserFunction< VectorType, double >      ValueFunctionType;
  typedef UserFunction< VectorType, VectorType >  DerivativeFunctionType;

  typedef UserFunction< double, double >    OptimizerValueFunctionType;
  typedef UserFunction< double, double >    OptimizerDerivativeFunctionType;

  /** Return the type of this object. */
  tubeTypeMacro( OptimizerND );

  /** Constructor. */
  OptimizerND( void );

  /** Constructor. */
  OptimizerND( unsigned int dimension,
    ValueFunctionType::Pointer funcValND,
    DerivativeFunctionType::Pointer funcDerivND,
    Optimizer1D::Pointer optimizer1D );

  /** Destructor. */
  virtual ~OptimizerND( void );

  tubeGetMacro( MaxIterations, unsigned int );

  virtual void SetMaxIterations( unsigned int maxIterations );

  tubeGetMacro( MaxLineSearches, unsigned int );

  tubeSetMacro( MaxLineSearches, unsigned int );

  tubeGetMacro( SearchForMin, bool );

  virtual void SetSearchForMin( bool searchForMin );

  tubeBooleanMacro( SearchForMin );

  tubeGetMacro( Tolerance, double );

  virtual void SetTolerance( double tolerance );

  tubeGetMacro( XMax, VectorType );

  tubeSetMacro( XMax, VectorType );

  tubeGetMacro( XMin, VectorType );

  tubeSetMacro( XMin, VectorType );

  tubeGetMacro( XStep, VectorType );

  tubeSetMacro( XStep, VectorType );

  bool Extreme( VectorType & x, double * xVal );

  bool Extreme( VectorType & x, double * xVal, unsigned int n,
    MatrixType & directions );

  double FuncDeriv( double x );

  double FuncVal( double x );

  void Use( unsigned int dimension,
    ValueFunctionType::Pointer funcValND,
    DerivativeFunctionType::Pointer funcDerivND,
    Optimizer1D::Pointer optimizer1D );

protected:

  /** Print out information about this object. */
  void PrintSelf( std::ostream & os, Indent indent ) const override;

  unsigned int                              m_Dimension;
  VectorType                                m_XMin;
  VectorType                                m_XMax;
  VectorType                                m_XStep;
  VectorType                                m_X0;
  VectorType                                m_X0Dir;
  VectorType                                m_X0Temp;
  bool                                      m_SearchForMin;
  double                                    m_Tolerance;
  unsigned int                              m_MaxIterations;
  unsigned int                              m_MaxLineSearches;
  OptimizerValueFunctionType::Pointer       m_Optimizer1DVal;
  OptimizerDerivativeFunctionType::Pointer  m_Optimizer1DDeriv;
  Optimizer1D::Pointer                      m_Optimizer1D;
  ValueFunctionType::Pointer                m_FuncValND;
  DerivativeFunctionType::Pointer           m_FuncDerivND;

private:

  // Copy constructor not implemented.
  OptimizerND( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

}; // End class OptimizerND

} // End namespace tube

#endif // End !defined( __tubeOptimizerND_h )
