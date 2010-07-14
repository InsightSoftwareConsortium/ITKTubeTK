/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkOptimizerND.cxx,v $
  Language:  C++
  Date:      $Date: 2003/12/16 14:07:06 $
  Version:   $Revision: 1.5 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkOptimizerND.h"
#include "itkOptimizer1D.h"
#include <iostream>
#include <cmath>

namespace itk
{


class OptValFuncND : public UserFunc<double, double> {
    OptimizerND * cOpt;
public:
  OptValFuncND(OptimizerND * newOpt)
  {
     cOpt = newOpt;
  };

  double value(double x)
  {
     return cOpt->funcVal(x);
  };
};

class OptDerivFuncND : public UserFunc<double, double> {
    OptimizerND * cOpt;
public:
  OptDerivFuncND(OptimizerND * newOpt)
  {
    cOpt = newOpt;
  };
  double value(double x)
  {
    return cOpt->funcDeriv(x);
  }
};


OptimizerND::OptimizerND()
{
  cDefined = false;
  cOpt1DVal = new OptValFuncND(this);
  cOpt1DDeriv = new OptDerivFuncND(this);
}


OptimizerND::OptimizerND(int newNDims, UserFunc<VectorType *, double> * newFuncValND, UserFunc<VectorType *, VectorType &> * newFuncDerivND, Optimizer1D *newOpt1D)
{
  cDefined = false;
  cOpt1DVal = new OptValFuncND(this);
  cOpt1DDeriv = new OptDerivFuncND(this);
  use(newNDims, newFuncValND, newFuncDerivND, newOpt1D);
}

OptimizerND::~OptimizerND()
{
  if(cDefined) 
  {
    delete cXMin;
    delete cXMax;
    delete cXStep;
    delete cX0;
    delete cX0Dir;
    delete cX0Temp;
  }
  delete cOpt1DVal;
  delete cOpt1DDeriv;
}

void OptimizerND::use(int newNDims, UserFunc<VectorType*, double> * newFuncValND, UserFunc<VectorType*, VectorType &> * newFuncDerivND, Optimizer1D *newOpt1D)
{
  if(cDefined) 
  {
    delete cXMin;
    delete cXMax;
    delete cXStep;
    delete cX0;
    delete cX0Dir;
    delete cX0Temp;
  }

  cNDims = newNDims;
  cXMin = new VectorType(cNDims, 0.0);
  cXMax = new VectorType(cNDims, 1.0);
  cXStep = new VectorType(cNDims, 0.01);
  cX0 = new VectorType(cNDims, 0.0);
  cX0Dir = new VectorType(cNDims, 1/sqrt((float)cNDims));
  cX0Temp = new VectorType(cNDims, 0.0);

  cFuncValND = newFuncValND;
  cFuncDerivND = newFuncDerivND;

  cSearchForMin = true;
  cTolerance = 0.0001;
  cMaxIterations = 300;

  cOpt1D = newOpt1D;
  cOpt1D->use(cOpt1DVal, cOpt1DDeriv);
  cDefined = true;
}


OptimizerND::VectorType & OptimizerND::xMin(void)
{
  return *cXMin;
}

void OptimizerND::xMin(VectorType & newXMinn)
{
   (*cXMin) = newXMinn;
}

OptimizerND::VectorType & OptimizerND::xMax(void)
{
  return *cXMax;
}

void OptimizerND::xMax(VectorType & newXMaxx)
{
  (*cXMax) = newXMaxx;
}

OptimizerND::VectorType & OptimizerND::xStep(void)
{
  return *cXStep;
}

void OptimizerND::xStep(VectorType & newXStepp)
{
  (*cXStep) = newXStepp;
}

double OptimizerND::tolerance(void)
{
  return cTolerance;
}


void OptimizerND::tolerance(double newTolerance)
{
  cOpt1D->tolerance(newTolerance);
  cTolerance = newTolerance;
}

unsigned int OptimizerND::maxIterations(void)
{
  return cMaxIterations;
}

void OptimizerND::maxIterations(unsigned int newMaxIterations)
{
  cOpt1D->maxIterations(newMaxIterations);
  cMaxIterations = newMaxIterations;
}


bool OptimizerND::searchForMin(void)
{
  return cSearchForMin;
}


void OptimizerND::searchForMin(bool newSearchForMin)
{
  cOpt1D->searchForMin(newSearchForMin);
  cSearchForMin = newSearchForMin;
}


double OptimizerND::funcVal(double a)
{
  *cX0Temp=ComputeLineStep(*cX0, a, *cX0Dir);
  return cFuncValND->value(cX0Temp);
}


double OptimizerND::funcDeriv(double a)
{
  *cX0Temp =ComputeLineStep(*cX0, a, *cX0Dir);
  return dot_product(cFuncDerivND->value(cX0Temp),*cX0Dir);
}


bool OptimizerND::extreme(VectorType & x, double * xVal)
{
  (*cX0) = x;
  double a = 1;
  double xmin, xmax, xstep;
  unsigned int count = 0;
  while(fabs(a)>cTolerance) 
  {
    (*cX0Dir) = cFuncDerivND->value(cX0);
    if(cX0Dir->magnitude()<cTolerance*cTolerance) 
    {
      a = 0;
      break;
    }
    cX0Dir->normalize();
    //Normalize(*cX0Dir);
    xmin = dot_product((*cXMin) - (*cX0), *cX0Dir);
    xmax = dot_product((*cXMax) - (*cX0), *cX0Dir);
    
    xstep = fabs(dot_product(*cXStep, *cX0Dir));
    if(xmax<xmin) 
    {
      a = xmin;
      xmin = xmax;
      xmax = a;
    }
    a = 0;
    cOpt1D->xMin(xmin);
    cOpt1D->xMax(xmax);
    cOpt1D->xStep(xstep);

    cOpt1D->extreme(&a, xVal);
    *cX0 = ComputeLineStep(*cX0, a, *cX0Dir);
    if(count++>2*cNDims)
    {
      break;
    }
  }
  x = (*cX0);
  if(fabs(a)>cTolerance) 
  {
    std::cout << "Scanned " << count << " directions without convergence - aborting" << std::endl;
    return false;
  }

  return true;
}

    

bool OptimizerND::extreme(VectorType & x, double * xVal, unsigned int n, MatrixType &dirs)
{
  (*cX0) = x;
  double a;
  double xmin, xmax, xstep;
  unsigned int i, j;
  for(i=0; i<n; i++) 
  {
    for(j=0; j<x.size(); j++)
    {
      (*cX0Dir)(j) = dirs.get(j,i);
    }

    //Normalize(*cX0Dir);
    cX0Dir->normalize();

    xmin = 0;
    xmax = 0;
    
    for(j=0; j<x.size(); j++)
    {
      if((*cX0Dir)(j)>0)
      {
        if((*cXMax)(j)-(*cX0)(j)>(*cXMin)(j)-(*cX0)(j))
        {
          xmin += (*cX0Dir)(j) * ((*cXMin)(j)-(*cX0)(j));
          xmax += (*cX0Dir)(j) * ((*cXMax)(j)-(*cX0)(j));
        }
        else
        {
          xmin += (*cX0Dir)(j) * ((*cXMax)(j)-(*cX0)(j));
          xmax += (*cX0Dir)(j) * ((*cXMin)(j)-(*cX0)(j));
        }
      }
      else
      {
        if((*cXMax)(j)-(*cX0)(j)<(*cXMin)(j)-(*cX0)(j))
        {
          xmin += (*cX0Dir)(j) * ((*cXMin)(j)-(*cX0)(j));
          xmax += (*cX0Dir)(j) * ((*cXMax)(j)-(*cX0)(j));
        }
        else
        {
          xmin += (*cX0Dir)(j) * ((*cXMax)(j)-(*cX0)(j));
          xmax += (*cX0Dir)(j) * ((*cXMin)(j)-(*cX0)(j));
        }
      }
    }

    xstep=0;
    for(j=0; j<x.size(); j++)
    {
      if((*cX0Dir)(j)>0)
      {        
        xstep += (*cXStep)(j) * (*cX0Dir)(j);
      }
      else
      {
        xstep -= (*cXStep)(j) * (*cX0Dir)(j);
      }
    }

    if(xmax<xmin) 
    {
      a = xmin;
      xmin = xmax;
      xmax = a;
    }
    a = 0;
    cOpt1D->xMin(xmin);
    cOpt1D->xMax(xmax);
    cOpt1D->xStep(xstep);
    cOpt1D->extreme(&a, xVal);
    *cX0=ComputeLineStep(*cX0, a, *cX0Dir);

  }
  x = (*cX0);
  return true;
}

}; // namespace ITK

