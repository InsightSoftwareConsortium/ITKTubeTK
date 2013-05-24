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

#include "TubeRegistrator.h"

#include <ctime>

#include <vnl/algo/vnl_matrix_inverse.h>

using namespace VTREE;
using namespace VHTK;

TubeRegistratorPoint::
TubeRegistratorPoint( void )
:m_X(3), m_V1(3), m_m_V2(3),
 m_XT(3), m_V1T(3), m_V2T(3), m_DXT(3)
{
  m_VAL = 0;
  m_X = 0;
  m_R = 0;
  m_V1 = 0;
  m_V2 = 0;
  m_RN = 0;
  m_MN = 0;
  m_W = 0;
  m_XT = 0;
  m_V1T = 0;
  m_V2T = 0;
  m_DXT = 0;
}

TubeRegistrator::
TubeRegistrator( void )
: m_BiasV(3,3), m_BiasVI(3,3)
{
  m_Kappa = 1;
  m_TubeNet = NULL;
  m_RegPoints.clear();
  m_Im = NULL;
  m_ImMin = 0;
  m_ImRange = 0;
  m_Sampling = 100;
  m_Count = 0;
  m_Weight = 0;
  m_BiasV = 0;
  m_BiasVI = 0;
  m_Metric = 0;
  m_RegImThresh = 0;
}

void TubeRegistrator::
SetTubeNet(TubeNet * tubes)
{
  m_TubeNet = tubes;
}

void TubeRegistrator::
SetImage(Image3D<short> * im)
{
  m_Im = im;
  m_ImOp.use(m_Im);

  m_ImMin = m_Im->dataMin();
  m_ImRange = m_Im->dataMax()-m_ImMin;
}

void TubeRegistrator::
SetSampling(int newSampling)
{
  m_Sampling = newSampling;
}

int TubeRegistrator::
GetSampling( void )
{
  return m_Sampling;
}

int TubeRegistrator::
SetNumSamples( void )
{
  return m_Count;
}

std::list<TubeRegistratorPoint *> * TubeRegistrator::
GetSamples( void )
{
  return & m_RegPoints;
}

void TubeRegistrator::
SetKappa(double kappa)
{
  m_Kappa = kappa;
}

void TubeRegistrator::
SetImThresh(double newRegThresh)
{
  m_RegImThresh = newRegThresh;
}

double TubeRegistrator::
GetImThresh( void )
{
  return m_RegImThresh;
}

void TubeRegistrator::
MetricPreProc( void )
{
  if(m_Im == NULL)
    {
    return;
    }

  if(m_TubeNet == NULL)
    {
    return;
    }

  std::list<Tube *>::iterator i;
  std::list<TubePoint *>::iterator j;
  std::list<TubeRegistratorPoint *>::iterator l;
  vnl_vector<double> tV(3);
  vnl_matrix<double> tM(3,3);

  for(l=m_RegPoints.begin(); l!=m_RegPoints.end(); ++l)
    {
    delete * l;
    }
  m_RegPoints.clear();

  m_Count = 0;
  m_Weight = 0;
  m_BiasV = 0;

  vnl_vector<double> tv(3);
  double cX = 0;
  double cY = 0;
  double cZ = 0;
  TubeRegistratorPoint * tPnt;
  for(i=m_TubeNet->tubes()->begin(); i!=m_TubeNet->tubes()->end(); ++i)
    {
    int tubeSize = (*i)->points()->size();
    if(tubeSize>m_Sampling)
      {
      int skipped = 0;
      for(j=(*i)->points()->begin(); j!=(*i)->points()->end(); ++j)
        {
        while(skipped++%(m_Sampling/2) != 0 && j!=(*i)->points()->end())
          {
          ++j;
          }
        if(j!=(*i)->points()->end() && skipped+10<tubeSize
           && (*j)->ridgeness()>0.2
           && (*j)->medialness()>0.02)
          {
          tPnt = new TubeRegistratorPoint;
          tPnt->m_X(0) = (*((*j)->m_X()))(1);
          tPnt->m_X(1) = (*((*j)->m_X()))(2);
          tPnt->m_X(2) = (*((*j)->m_X()))(3);
          tPnt->m_R = (*j)->m_R();
          tPnt->m_RN = (*j)->ridgeness();
          tPnt->m_MN = (*j)->medialness();
          tPnt->m_V1(0) = (*((*j)->m_V1()))(1);
          tPnt->m_V1(1) = (*((*j)->m_V1()))(2);
          tPnt->m_V1(2) = (*((*j)->m_V1()))(3);
          tPnt->m_V2(0) = (*((*j)->m_V2()))(1);
          tPnt->m_V2(1) = (*((*j)->m_V2()))(2);
          tPnt->m_V2(2) = (*((*j)->m_V2()))(3);

          tPnt->m_W = 2/(1+vcl_exp(-4*tPnt->m_R))-1;

          m_RegPoints.push_back(tPnt);

          cX += tPnt->m_W*tPnt->m_X(0);
          cY += tPnt->m_W*tPnt->m_X(1);
          cZ += tPnt->m_W*tPnt->m_X(2);

          tM = outer_product(tPnt->m_V1T, tPnt->m_V1T);
          tM = tM + outer_product(tPnt->m_V2T, tPnt->m_V2T);
          tM = tPnt->m_W * tM;
          m_BiasV = m_BiasV + tM;

          m_Weight += tPnt->m_W;
          m_Count++;
          }
        while(skipped++%(m_Sampling/2) != 0 && j!=(*i)->points()->end())
          {
          ++j;
          }
        if(j==(*i)->points()->end())
          {
          --j;
          }
        }
      }
    }
  m_BiasV = (1.0/m_Weight) * m_BiasV;

  SetCenter(cX/m_Weight, cY/m_Weight, cZ/m_Weight);

  std::cout << "Sampling = " << m_Sampling << std::endl;
  std::cout << "Num Samples = " << m_Count << std::endl;
  std::cout << "Center = " << cX/m_Weight << ", " << cY/m_Weight << ", " << cZ/m_Weight << std::endl;
}

double TubeRegistrator::
Metric( void )
{
  long c0 = std::clock();

  std::list<TubeRegistratorPoint *>::iterator j;
  static TNT::Vector<double> xTV(3);
  static TNT::Vector<double> dXTV(3);

  bool recalc = Updated();
  if(!recalc)
    {
    std::cout << "Metric: !recalc" << std::endl;
    std::cout << "Time = " << (std::clock()-c0)/(double)CLOCKS_PER_SEC << std::endl;
    return m_Metric;
    }

  m_Metric = 0;
  m_Count = 0;
  m_Weight = 0;

  if(m_Kappa == 0)
    {
    for(j=m_RegPoints.begin(); j!=m_RegPoints.end(); ++j)
      {
      TransformPoint(&((*j)->m_X), &((*j)->m_XT));
      if((*j)->m_XT(0)>=0 && (*j)->m_XT(0)<m_Im->dimSizeX()-1 &&
         (*j)->m_XT(1)>=0 && (*j)->m_XT(1)<m_Im->dimSizeY()-1 &&
         (*j)->m_XT(2)>=0 && (*j)->m_XT(2)<m_Im->dimSizeZ()-1)
        {
        m_Weight += 1;
        m_Count++;
        (*j)->m_VAL = (*m_Im)((int)(*j)->m_XT(0),
                              (int)(*j)->m_XT(1),
                              (int)(*j)->m_XT(2));
        m_Metric += (*j)->m_VAL;
        }
      }
    }
  else
    {
    for(j=m_RegPoints.begin(); j!=m_RegPoints.end(); ++j)
      {
      TransformPoint(&((*j)->m_X), &((*j)->m_XT));
      if((*j)->m_XT(0)>0 && (*j)->m_XT(0)<m_Im->dimSizeX()-2 &&
         (*j)->m_XT(1)>0 && (*j)->m_XT(1)<m_Im->dimSizeY()-2 &&
         (*j)->m_XT(2)>0 && (*j)->m_XT(2)<m_Im->dimSizeZ()-2 &&
         (*m_Im)((*j)->m_XT(0), (*j)->m_XT(1), (*j)->m_XT(2)) > m_RegImThresh)
        {
        xTV(1) = (*j)->m_XT(0);
        xTV(2) = (*j)->m_XT(1);
        xTV(3) = (*j)->m_XT(2);
        m_Weight += (*j)->m_W;
        m_Count++;
        double opR = (*j)->r;
        if(opR<0.5)
          opR = 0.5;
        m_ImOp.opScale(opR * m_Kappa);
        //(*j)->m_VAL = (*j)->m_W * m_ImOp.valueIso(xTV);
        (*j)->m_VAL = (*j)->m_W * m_ImOp.valueIsoD(xTV, dXTV);
        (*j)->m_DXT(0) = dXTV(1);
        (*j)->m_DXT(1) = dXTV(2);
        (*j)->m_DXT(2) = dXTV(3);
        m_Metric += (*j)->m_VAL;
        }
      }
    }

  if(m_Weight == 0 || m_ImRange == 0)
    {
    m_Metric = 10;
    }
  else
    {
    m_Metric = (m_Metric/m_Weight-m_ImMin)/m_ImRange;
    }

  //std::cout << "Time = " << (std::clock()-c0)/(double)CLOCKS_PER_SEC << std::endl;

  return m_Metric;
}

double TubeRegistrator::
Metrim_Deriv(double *dX, double *dY, double *dZ,
             double *dA, double *dB, double *dG)
{

  long c0 = std::clock();

  double opR;

  bool recalc = Updated();

  if(recalc)
    {
    m_Metric = 0;
    m_Weight = 0;
    m_Count = 0;
    }
  else
    std::cout << "Metrim_Deriv: !recalc" << std::endl;

  *dX = 0;
  *dY = 0;
  *dZ = 0;
  *dA = 0;
  *dB = 0;
  *dG = 0;

  m_BiasV = 0;

  double dXProj1;
  double dXProj2;

  std::list<TubeRegistratorPoint *>::iterator j;
  TNT::Vector<double> xTV(3);
  TNT::Vector<double> dXTV(3);
  vnl_vector<double> tV(3);
  vnl_matrix<double> tM(3,3);

  double tDA;
  double tDB;
  double tDG;

  //FILE * fp = std::fopen("test.dat", "m_W");
  for(j=m_RegPoints.begin(); j!=m_RegPoints.end(); ++j)
    {
    if(recalc)
      TransformPoint(&((*j)->m_X), &((*j)->m_XT));
    if((*j)->m_XT(0)>0 && (*j)->m_XT(0)<m_Im->dimSizeX()-2 &&
       (*j)->m_XT(1)>0 && (*j)->m_XT(1)<m_Im->dimSizeY()-2 &&
       (*j)->m_XT(2)>0 && (*j)->m_XT(2)<m_Im->dimSizeZ()-2 &&
       (*m_Im)((*j)->m_XT(0), (*j)->m_XT(1), (*j)->m_XT(2)) > m_RegImThresh)
      {
      if(recalc)
        {
        m_Weight += (*j)->m_W;
        m_Count++;
        opR = (*j)->r;
        if(opR<0.5)
          opR = 0.5;
        m_ImOp.opScale(opR * m_Kappa);
        xTV(1) = (*j)->m_XT(0);
        xTV(2) = (*j)->m_XT(1);
        xTV(3) = (*j)->m_XT(2);
        (*j)->m_VAL = (*j)->m_W * m_ImOp.valueIsoD(xTV, dXTV);
        (*j)->m_DXT(0) = dXTV(1);
        (*j)->m_DXT(1) = dXTV(2);
        (*j)->m_DXT(2) = dXTV(3);
        m_Metric += (*j)->m_VAL;
        }

      TransformCoVector(&((*j)->m_V1), &((*j)->m_V1T));
      TransformCoVector(&((*j)->m_V2), &((*j)->m_V2T));

      dXProj1 = dot_product((*j)->m_DXT, (*j)->m_V1T);
      dXProj2 = dot_product((*j)->m_DXT, (*j)->m_V2T);
      //std::fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
      //        (*j)->m_X(0), (*j)->m_X(1), (*j)->m_X(2),
      //        (*j)->m_XT(0), (*j)->m_XT(1), (*j)->m_XT(2),
      //        (*j)->m_VAL,
      //        (*j)->m_DXT(0), (*j)->m_DXT(1), (*j)->m_DXT(2),
      //        (*j)->m_V1T(0), (*j)->m_V1T(1), (*j)->m_V1T(2),
      //        (*j)->m_V2T(0), (*j)->m_V2T(1), (*j)->m_V2T(2),
      //        dXProj1, dXProj2, (*j)->m_R, (*j)->m_RN, (*j)->m_MN);
      (*j)->m_DXT = (dXProj1*(*j)->m_V1T + dXProj2*(*j)->m_V2T);

      tM = outer_product((*j)->m_V1T, (*j)->m_V1T);
      tM = tM + outer_product((*j)->m_V2T, (*j)->m_V2T);
      tM = (*j)->m_W * tM;
      m_BiasV = m_BiasV + tM;

      *dX += (*j)->m_W * (*j)->m_DXT(0);
      *dY += (*j)->m_W * (*j)->m_DXT(1);
      *dZ += (*j)->m_W * (*j)->m_DXT(2);
      }
    }
  //std::fclose(fp);

  m_BiasVI = vnl_matrix_inverse<double>(m_BiasV).inverse();
  tV(0) = *dX;
  tV(1) = *dY;
  tV(2) = *dZ;
  tV = tV * m_BiasVI;
  *dX = tV(0);
  *dY = tV(1);
  *dZ = tV(2);

  if(m_Weight == 0)
    {
    m_BiasV = 0;
    *dA = 0;
    *dB = 0;
    *dG = 0;
    m_Metric = 10;
    return 10;
    }
  else
    {
    m_BiasV = 1.0/m_Weight * m_BiasV;
    }
  m_BiasVI = vnl_matrix_inverse<double>(m_BiasV).inverse();

  for(j=m_RegPoints.begin(); j!=m_RegPoints.end(); ++j)
    {
    if((*j)->m_XT(0)>=0 && (*j)->m_XT(0)<m_Im->dimSizeX()-1 &&
       (*j)->m_XT(1)>=0 && (*j)->m_XT(1)<m_Im->dimSizeY()-1 &&
       (*j)->m_XT(2)>=0 && (*j)->m_XT(2)<m_Im->dimSizeZ()-1)
      {
      (*j)->m_DXT = (*j)->m_DXT * m_BiasVI;
      //(*j)->m_DXT(0) -= (*dX);
      //(*j)->m_DXT(1) -= (*dY);
      //(*j)->m_DXT(2) -= (*dZ);
      GetDeltaAngles(&((*j)->m_XT), &((*j)->m_DXT), &tDA, &tDB, &tDG);
      //std::cout << "m_X = " << (*j)->m_XT(0) << ", " << (*j)->m_XT(1) << ", " << (*j)->m_XT(2);
      //std::cout << " : dA = " << tDA << ", " << tDB << ", " << tDG << std::endl;
      *dA += (*j)->m_W * tDA;
      *dB += (*j)->m_W * tDB;
      *dG += (*j)->m_W * tDG;
      }
    }
  *dA /= m_Weight;
  *dB /= m_Weight;
  *dG /= m_Weight;

  m_Metric = (m_Metric/m_Weight-m_ImMin)/m_ImRange;

  std::cout << "Time = " << (std::clock()-c0)/(double)CLOCKS_PER_SEC << std::endl;

  return m_Metric;
}

//
//
//

#include <MathLib/OptimizerND.h>
#include <MathLib/OptParabolicFit1D.h>
#include <MathLib/OptBrent1D.h>
#include <MathLib/OptGrad1D.h>
#include <MathLib/UserFunc.h>

double offsetUnit=2.5;
double rotUnit=0.1;
double offsetDUnit=3.5;
double rotDUnit=0.1;

class mVMetricCost : public UserFunc<TNT::Vector<double> *, double>
{
public:
  mVMetricCost( void )
    : m_RegOp(0)
    {
    }
  void use(TubeRegistrator * newRegOp)
    {
    m_RegOp = newRegOp;
    }
  double value(TNT::Vector<double> * m_X)
    {
    m_RegOp->SetOffset((*m_X)(1)*offsetUnit, (*m_X)(2)*offsetUnit, (*m_X)(3)*offsetUnit);
    m_RegOp->SetAngles((*m_X)(4)*rotUnit, (*m_X)(5)*rotUnit, (*m_X)(6)*rotUnit);
    /*
    std::cout << "V "
              << (*m_X)(1) << " (" << (*m_X)(1)*offsetUnit << ")  "
              << (*m_X)(2) << " (" << (*m_X)(2)*offsetUnit << ")  "
              << (*m_X)(3) << " (" << (*m_X)(3)*offsetUnit << ")  " << std::endl;
    std::cout << "  "
              << (*m_X)(4) << " (" << (*m_X)(4)*rotUnit << ")  "
              << (*m_X)(5) << " (" << (*m_X)(5)*rotUnit << ")  "
              << (*m_X)(6) << " (" << (*m_X)(6)*rotUnit << ")  " << std::endl;
    */
    double tf = m_RegOp->Metric();
    //std::cout << "   v = " << tf << std::endl;
    return tf;
    }

protected:
  TubeRegistrator  * m_RegOp;

}; // End class mVMetricCost

class mVMetricDeriv : public UserFunc<TNT::Vector<double> *, TNT::Vector<double> *>
{
public:
  mVMetricDeriv( void )
    : m_RegOp(0), m_D(6, 0.0)
    {
    }
  void use(TubeRegistrator * newRegOp)
    {
    m_RegOp = newRegOp;
    }
  TNT::Vector<double> * value(TNT::Vector<double> * m_X)
    {
    m_RegOp->SetOffset((*m_X)(1)*offsetUnit, (*m_X)(2)*offsetUnit, (*m_X)(3)*offsetUnit);
    m_RegOp->SetAngles((*m_X)(4)*rotUnit, (*m_X)(5)*rotUnit, (*m_X)(6)*rotUnit);
    /*
    std::cout << "D "
             << (*m_X)(1) << " (" << (*m_X)(1)*offsetUnit << ")  "
             << (*m_X)(2) << " (" << (*m_X)(2)*offsetUnit << ")  "
             << (*m_X)(3) << " (" << (*m_X)(3)*offsetUnit << ")  " << std::endl;
    std::cout << "  "
             << (*m_X)(4) << " (" << (*m_X)(4)*rotUnit << ")  "
             << (*m_X)(5) << " (" << (*m_X)(5)*rotUnit << ")  "
             << (*m_X)(6) << " (" << (*m_X)(6)*rotUnit << ")  " << std::endl;
    */
    m_RegOp->Metrim_Deriv(&m_D(1), &m_D(2), &m_D(3),
                        &m_D(4), &m_D(5), &m_D(6));
    /*
    std::cout << "  "
             << m_D(1) << " (" << m_D(1)/offsetDUnit << ")  "
             << m_D(2) << " (" << m_D(2)/offsetDUnit << ")  "
             << m_D(3) << " (" << m_D(3)/offsetDUnit << ")  " << std::endl;
    std::cout << "  "
             << m_D(4) << " (" << m_D(4)/rotDUnit << ")  "
             << m_D(5) << " (" << m_D(5)/rotDUnit << ")  "
             << m_D(6) << " (" << m_D(6)/rotDUnit << ")  " << std::endl;
    */
    m_D(1) /= offsetDUnit;
    m_D(2) /= offsetDUnit;
    m_D(3) /= offsetDUnit;
    m_D(4) /= rotDUnit;
    m_D(5) /= rotDUnit;
    m_D(6) /= rotDUnit;
    return &m_D;
    }

protected:
  TNT::Vector<double>  m_D;
  TubeRegistrator    * m_RegOp;

}; // End class mVMetricDeriv

bool TubeRegistrator::Fit( void )
{
  mVMetricCost * func = new mVMetricCost;
  func->use(this);

  mVMetricDeriv * funcD = new mVMetricDeriv;
  funcD->use(this);

  OptParabolicFit1D op1D;
  OptimizerND op(6, func, funcD, &op1D);
  op.searchForMin(false);

  double x = (*mO)(0);
  double y = (*mO)(1);
  double z = (*mO)(2);
  double a = mAlpha;
  double b = mBeta;
  double g = mGamma;

  TNT::Vector<double> m(6);
  m = -100;
  op.xMin(&m);
  m = 100;
  op.xMax(&m);
  m = 1.0;
  op.xStep(&m);

  std::cout << "Initial State = " << x << ", " << y << ", " << z << std::endl;
  std::cout << "                " << a << ", " << b << ", " << g << std::endl;
  std::cout << "   m_VAL = " << func->value(&m) << std::endl;

  double xVal;
  m(1) = x/offsetUnit;
  m(2) = y/offsetUnit;
  m(3) = z/offsetUnit;
  m(4) = a/rotUnit;
  m(5) = b/rotUnit;
  m(6) = g/rotUnit;

  m_Sampling *= 2;
  //m_Kappa *= 1.5;
  this->MetricPreProc();

  std::cout << "Bounding min..." << std::endl;
  op.bound(&m, &xVal);

  std::cout << "Min bound at..."
           << m(1) << " (" << m(1)*offsetUnit << ")  "
           << m(2) << " (" << m(2)*offsetUnit << ")  "
           << m(3) << " (" << m(3)*offsetUnit << ")  " << std::endl;
  std::cout << "  "
           << m(4) << " (" << m(4)*rotUnit << ")  "
           << m(5) << " (" << m(5)*rotUnit << ")  "
           << m(6) << " (" << m(6)*rotUnit << ")  " << std::endl;

  m_Sampling /= 2;
  //m_Kappa /= 1.5;
  MetricPreProc();
  xVal = func->value(&m);
  op.maxIterations(20);
  op.tolerance(0.0001);
  op.extreme(&m, &xVal);

  x = m(1)*offsetUnit;
  y = m(2)*offsetUnit;
  z = m(3)*offsetUnit;
  a = m(4)*rotUnit;
  b = m(5)*rotUnit;
  g = m(6)*rotUnit;

  std::cout << "Final State = " << x << ", " << y << ", " << z << std::endl;
  std::cout << "              " << a << ", " << b << ", " << g << std::endl;
  std::cout << "   m_VAL = " << xVal << std::endl;

  (*mO)(0) = x;
  (*mO)(1) = y;
  (*mO)(2) = z;
  mAlpha = a;
  mBeta = b;
  mGamma = g;

  return true;
}
