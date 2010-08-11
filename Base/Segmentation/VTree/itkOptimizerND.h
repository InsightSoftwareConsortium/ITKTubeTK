/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkOptimizerND.h,v $
  Language:  C++
  Date:      $Date: 2003/01/13 19:59:25 $
  Version:   $Revision: 1.4 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

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

  bool         searchForMin( void );
  void         searchForMin( bool newSearchForMin );

  double       funcVal( double x );
  double       funcDeriv( double x );

  bool         extreme( vnl_vector<double> &x, double *xVal );
  bool         extreme( vnl_vector<double> &x, double *xVal, unsigned int n,
                        MatrixType &dirs );

protected :
       
  bool         cDefined;
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

  UserFunc< double, double >             * cOpt1DVal;
  UserFunc< double, double >             * cOpt1DDeriv;

  Optimizer1D                            * cOpt1D;

  UserFunc< vnl_vector<double>, double >               * cFuncValND;
  UserFunc< vnl_vector<double>, vnl_vector<double> > * cFuncDerivND;
  
};



}; // end namespace itk



#endif

