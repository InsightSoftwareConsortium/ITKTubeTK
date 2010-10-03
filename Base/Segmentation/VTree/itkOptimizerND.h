/*=========================================================================

Library:   TubeTK/VTree3D

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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
#ifndef __itkOptimizerND_h
#define __itkOptimizerND_h

#include "itkUserFunc.h"
#include "itkOptimizer1D.h"

#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>

#include "itkMatrixMath.h"

namespace itk 
{


class OptimizerND
{

public :

  /** Typedef for the matrix type used */
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

protected :
       
  unsigned int cNDims;
  
  vnl_vector<double>   cXMin;
  vnl_vector<double>   cXMax;
  vnl_vector<double>   cXStep;
  vnl_vector<double>   cX0;
  vnl_vector<double>   cX0Dir;
  vnl_vector<double>   cX0Temp;

  bool         cSearchForMin;
  double       cTolerance;
  unsigned int cMaxIterations;
  unsigned int cMaxLineSearches;

  UserFunc< double, double >             * cOpt1DVal;
  UserFunc< double, double >             * cOpt1DDeriv;

  Optimizer1D                            * cOpt1D;

  UserFunc< vnl_vector<double>, double >               * cFuncValND;
  UserFunc< vnl_vector<double>, vnl_vector<double> > * cFuncDerivND;
  
};



}; // end namespace itk



#endif

