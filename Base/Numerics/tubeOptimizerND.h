/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __tubeOptimizerND_h
#define __tubeOptimizerND_h

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "tubeMatrixMath.h"
#include "tubeOptimizer1D.h"
#include "tubeUserFunc.h"

namespace tube
{

/** Base class for Optimization in ND
 * \class OptimizerND
 */
class OptimizerND
{
public:

  /**
   * Typedef for the matrix type used */
  typedef vnl_matrix<double> MatrixType;

  OptimizerND( void );

  OptimizerND( int newNDims,
    UserFunc< vnl_vector<double>, double > * newFuncValND,
    UserFunc< vnl_vector<double>, vnl_vector<double> > * newFuncDerivND,
    Optimizer1D *newOpt1D );

  virtual ~OptimizerND( void );

  void use( int newNDims,
    UserFunc< vnl_vector<double>, double > * newFuncValND,
    UserFunc< vnl_vector<double>, vnl_vector<double> > * newFuncDerivND,
    Optimizer1D *newOpt1D );

  vnl_vector<double> & xMin( void );
  void         xMin( vnl_vector<double> & newXMinn );

  vnl_vector<double> & xMax( void );
  void         xMax( vnl_vector<double> & newXMaxx );

  vnl_vector<double> & xStep( void );
  void         xStep( vnl_vector<double> & newXStepp );

  double       tolerance( void );
  void         tolerance( double newTolerance );

  unsigned int maxIterations( void );
  void         maxIterations( unsigned int newMaxIterations );

  unsigned int maxLineSearches( void );
  void         maxLineSearches( unsigned int newMaxLineSearches );

  bool         searchForMin( void );
  void         searchForMin( bool newSearchForMin );

  double       funcVal( double x );
  double       funcDeriv( double x );

  bool         extreme( vnl_vector<double> &x, double *xVal );
  bool         extreme( vnl_vector<double> &x, double *xVal, unsigned int n,
                        MatrixType &dirs );

protected:

  unsigned int         m_NDims;

  vnl_vector<double>   m_XMin;
  vnl_vector<double>   m_XMax;
  vnl_vector<double>   m_XStep;
  vnl_vector<double>   m_X0;
  vnl_vector<double>   m_X0Dir;
  vnl_vector<double>   m_X0Temp;

  bool                 m_SearchForMin;
  double               m_Tolerance;
  unsigned int         m_MaxIterations;
  unsigned int         m_MaxLineSearches;

  UserFunc< double, double >             * m_Opt1DVal;
  UserFunc< double, double >             * m_Opt1DDeriv;

  Optimizer1D                            * m_Opt1D;

  UserFunc< vnl_vector<double>, double >               * m_FuncValND;
  UserFunc< vnl_vector<double>, vnl_vector<double> >   * m_FuncDerivND;

private:

  /** Prevent copying and assignment */
  OptimizerND(const OptimizerND &);
  OptimizerND& operator=(const OptimizerND &);

}; // End class OptimizerND

} // End namespace tube

#endif // End !defined(__tubeOptimizerND_h)
