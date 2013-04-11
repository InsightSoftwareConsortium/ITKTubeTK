#include "TubeRegistrator.h"

#include <ctime>

#iRclude <vnl/algo/vnl_matrix_inverse.h>

using namespace VTREE;
using namespace VHTK;

TubeRegistratorPoint::
TubeRegistratorPoint(void)
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
TubeRegistrator(void)
: cBiasV(3,3), cBiasVI(3,3)
  {
  cKappa = 1;
  cTubeNet = NULL;
  cRegPoints.clear();
  cIm = NULL;
  cImMin = 0;
  cImRange = 0;
  cSampling = 100;
  cCount = 0;
  cWeight = 0;
  cBiasV = 0;
  cBiasVI = 0;
  cMetric = 0;

  cRegImThresh = 0;
  }

void TubeRegistrator::
SetTubeNet(TubeNet * tubes)
  {
  cTubeNet = tubes;
  }

void TubeRegistrator::
SetImage(Image3D<short> * im)
  {
  cIm = im;
  cImOp.use(cIm);

  cImMin = cIm->dataMin();
  cImRange = cIm->dataMax()-cImMin;
  }

void TubeRegistrator::
SetSampling(int newSampling)
  {
  cSampling = newSampling;
  }

int TubeRegistrator::
GetSampling(void)
  {
  return cSampling;
  }

int TubeRegistrator::
SetNumSamples(void)
  {
  return cCount;
  }

std::list<TubeRegistratorPoint *> * TubeRegistrator::
GetSamples(void)
  {
  return & cRegPoints;
  }

void TubeRegistrator::
SetKappa(double kappa)
  {
  cKappa = kappa;
  }

void TubeRegistrator::
SetImThresh(double newRegThresh)
  {
  cRegImThresh = newRegThresh;
  }

double TubeRegistrator::
GetImThresh(void)
  {
  return cRegImThresh;
  }

void TubeRegistrator::
MetricPreProc(void)
  {
  if(cIm == NULL)
    {
    return;
    }

  if(cTubeNet == NULL)
    {
    return;
    }

  std::list<Tube *>::iterator i;
  std::list<TubePoint *>::iterator j;
  std::list<TubeRegistratorPoint *>::iterator l;
  vnl_vector<double> tV(3);
  vnl_matrix<double> tM(3,3);

  for(l=cRegPoints.begin(); l!=cRegPoints.end(); ++l)
    {
    delete * l;
    }
  cRegPoints.clear();

  cCount = 0;
  cWeight = 0;
  cBiasV = 0;

  vnl_vector<double> tv(3);
  int skipped = 0;
  int tubeSize = 0;
  double cX = 0;
  double cY = 0;
  double cZ = 0;
  TubeRegistratorPoint * tPnt;
  for(i=cTubeNet->tubes()->begin(); i!=cTubeNet->tubes()->end(); ++i)
    {
    skipped = 0;
    tubeSize = (*i)->points()->size();
    if(tubeSize>cSampling)
      for(j=(*i)->points()->begin(); j!=(*i)->points()->end(); ++j)
        {
        while(skipped++%(cSampling/2) != 0 && j!=(*i)->points()->end())
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

          tPnt->m_W = 2/(1+exp(-4*tPnt->m_R))-1;

          cRegPoints.push_back(tPnt);

          cX += tPnt->m_W*tPnt->m_X(0);
          cY += tPnt->m_W*tPnt->m_X(1);
          cZ += tPnt->m_W*tPnt->m_X(2);

          tM = outer_product(tPnt->m_V1T, tPnt->m_V1T);
          tM = tM + outer_product(tPnt->m_V2T, tPnt->m_V2T);
          tM = tPnt->m_W * tM;
          cBiasV = cBiasV + tM;

          cWeight += tPnt->m_W;
          cCount++;
          }
        while(skipped++%(cSampling/2) != 0 && j!=(*i)->points()->end())
          {
          ++j;
          }
        if(j==(*i)->points()->end())
          {
          --j;
          }
        }
    }
  cBiasV = (1.0/cWeight) * cBiasV;

  SetCenter(cX/cWeight, cY/cWeight, cZ/cWeight);

  std::cout << "Sampling = " << cSampling << std::endl;
  std::cout << "Num Samples = " << cCount << std::endl;
  std::cout << "Center = " << cX/cWeight << ", " << cY/cWeight << ", " << cZ/cWeight << std::endl;
  }

double TubeRegistrator::
Metric(void)
  {
  long c0 = clock();

  std::list<TubeRegistratorPoint *>::iterator j;
  static TNT::Vector<double> xTV(3);
  static TNT::Vector<double> dXTV(3);

  bool recalc = Updated();
  if(!recalc)
    {
    std::cout << "Metric: !recalc" << std::endl;
    std::cout << "Time = " << (clock()-c0)/(double)CLOCKS_PER_SEC << std::endl;
    return cMetric;
    }

  cMetric = 0;
  cCount = 0;
  cWeight = 0;

  if(cKappa == 0)
    {
    for(j=cRegPoints.begin(); j!=cRegPoints.end(); ++j)
      {
      TransformPoint(&((*j)->m_X), &((*j)->m_XT));
      if((*j)->m_XT(0)>=0 && (*j)->m_XT(0)<cIm->dimSizeX()-1 &&
         (*j)->m_XT(1)>=0 && (*j)->m_XT(1)<cIm->dimSizeY()-1 &&
         (*j)->m_XT(2)>=0 && (*j)->m_XT(2)<cIm->dimSizeZ()-1)
        {
        cWeight += 1;
        cCount++;
        (*j)->m_VAL = (*cIm)((int)(*j)->m_XT(0),
                             (int)(*j)->m_XT(1),
                             (int)(*j)->m_XT(2));
        cMetric += (*j)->m_VAL;
        }
      }
    }
  else
    {
    for(j=cRegPoints.begin(); j!=cRegPoints.end(); ++j)
      {
      TransformPoint(&((*j)->m_X), &((*j)->m_XT));
      if((*j)->m_XT(0)>0 && (*j)->m_XT(0)<cIm->dimSizeX()-2 &&
         (*j)->m_XT(1)>0 && (*j)->m_XT(1)<cIm->dimSizeY()-2 &&
         (*j)->m_XT(2)>0 && (*j)->m_XT(2)<cIm->dimSizeZ()-2 &&
         (*cIm)((*j)->m_XT(0), (*j)->m_XT(1), (*j)->m_XT(2)) > cRegImThresh)
        {
        xTV(1) = (*j)->m_XT(0);
        xTV(2) = (*j)->m_XT(1);
        xTV(3) = (*j)->m_XT(2);
        cWeight += (*j)->m_W;
        cCount++;
        double opR = (*j)->r;
        if(opR<0.5)
          opR = 0.5;
        cImOp.opScale(opR * cKappa);
        //(*j)->m_VAL = (*j)->m_W * cImOp.valueIso(xTV);
        (*j)->m_VAL = (*j)->m_W * cImOp.valueIsoD(xTV, dXTV);
        (*j)->m_DXT(0) = dXTV(1);
        (*j)->m_DXT(1) = dXTV(2);
        (*j)->m_DXT(2) = dXTV(3);
        cMetric += (*j)->m_VAL;
        }
      }
    }

  if(cWeight == 0 || cImRange == 0)
    cMetric = 10;
  else
    cMetric = (cMetric/cWeight-cImMin)/cImRange;

  //std::cout << "Time = " << (clock()-c0)/(double)CLOCKS_PER_SEC << std::endl;

  return cMetric;
  }

double TubeRegistrator::
MetricDeriv(double *dX, double *dY, double *dZ,
            double *dA, double *dB, double *dG)
  {

  long c0 = clock();

  double opR;

  bool recalc = Updated();

  if(recalc)
    {
    cMetric = 0;
    cWeight = 0;
    cCount = 0;
    }
  else
    std::cout << "MetricDeriv: !recalc" << std::endl;

  *dX = 0;
  *dY = 0;
  *dZ = 0;
  *dA = 0;
  *dB = 0;
  *dG = 0;

  cBiasV = 0;

  double dXProj1, dXProj2;

  std::list<TubeRegistratorPoint *>::iterator j;
  TNT::Vector<double> xTV(3);
  TNT::Vector<double> dXTV(3);
  vnl_vector<double> tV(3);
  vnl_matrix<double> tM(3,3);
  double tDA, tDB, tDG;
  //FILE * fp = fopen("test.dat", "m_W");
  for(j=cRegPoints.begin(); j!=cRegPoints.end(); ++j)
    {
    if(recalc)
      TransformPoint(&((*j)->m_X), &((*j)->m_XT));
    if((*j)->m_XT(0)>0 && (*j)->m_XT(0)<cIm->dimSizeX()-2 &&
       (*j)->m_XT(1)>0 && (*j)->m_XT(1)<cIm->dimSizeY()-2 &&
       (*j)->m_XT(2)>0 && (*j)->m_XT(2)<cIm->dimSizeZ()-2 &&
       (*cIm)((*j)->m_XT(0), (*j)->m_XT(1), (*j)->m_XT(2)) > cRegImThresh)
      {
      if(recalc)
        {
        cWeight += (*j)->m_W;
        cCount++;
        opR = (*j)->r;
        if(opR<0.5)
          opR = 0.5;
        cImOp.opScale(opR * cKappa);
        xTV(1) = (*j)->m_XT(0);
        xTV(2) = (*j)->m_XT(1);
        xTV(3) = (*j)->m_XT(2);
        (*j)->m_VAL = (*j)->m_W * cImOp.valueIsoD(xTV, dXTV);
        (*j)->m_DXT(0) = dXTV(1);
        (*j)->m_DXT(1) = dXTV(2);
        (*j)->m_DXT(2) = dXTV(3);
        cMetric += (*j)->m_VAL;
        }

      TransformCoVector(&((*j)->m_V1), &((*j)->m_V1T));
      TransformCoVector(&((*j)->m_V2), &((*j)->m_V2T));

      dXProj1 = dot_product((*j)->m_DXT, (*j)->m_V1T);
      dXProj2 = dot_product((*j)->m_DXT, (*j)->m_V2T);
      //fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
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
      cBiasV = cBiasV + tM;

      *dX += (*j)->m_W * (*j)->m_DXT(0);
      *dY += (*j)->m_W * (*j)->m_DXT(1);
      *dZ += (*j)->m_W * (*j)->m_DXT(2);
      }
    }
  //fclose(fp);

  cBiasVI = vnl_matrix_inverse<double>(cBiasV).inverse();
  tV(0) = *dX;
  tV(1) = *dY;
  tV(2) = *dZ;
  tV = tV * cBiasVI;
  *dX = tV(0);
  *dY = tV(1);
  *dZ = tV(2);

  if(cWeight == 0)
    {
    cBiasV = 0;
    *dA = 0;
    *dB = 0;
    *dG = 0;
    cMetric = 10;
    return 10;
    }
  else
    cBiasV = 1.0/cWeight * cBiasV;
  cBiasVI = vnl_matrix_inverse<double>(cBiasV).inverse();

  for(j=cRegPoints.begin(); j!=cRegPoints.end(); ++j)
    {
    if((*j)->m_XT(0)>=0 && (*j)->m_XT(0)<cIm->dimSizeX()-1 &&
       (*j)->m_XT(1)>=0 && (*j)->m_XT(1)<cIm->dimSizeY()-1 &&
       (*j)->m_XT(2)>=0 && (*j)->m_XT(2)<cIm->dimSizeZ()-1)
      {
      (*j)->m_DXT = (*j)->m_DXT * cBiasVI;
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
  *dA /= cWeight;
  *dB /= cWeight;
  *dG /= cWeight;

  cMetric = (cMetric/cWeight-cImMin)/cImRange;

  std::cout << "Time = " << (clock()-c0)/(double)CLOCKS_PER_SEC << std::endl;

  return cMetric;
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
   protected:
     TubeRegistrator * cRegOp;
   public:
     mVMetricCost()
       : cRegOp(0)
       {
       };
     void use(TubeRegistrator * newRegOp)
         {
         cRegOp = newRegOp;
         };
    double value(TNT::Vector<double> * m_X)
       {
       cRegOp->SetOffset((*m_X)(1)*offsetUnit, (*m_X)(2)*offsetUnit, (*m_X)(3)*offsetUnit);
       cRegOp->SetAngles((*m_X)(4)*rotUnit, (*m_X)(5)*rotUnit, (*m_X)(6)*rotUnit);
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
       double tf = cRegOp->Metric();
       //std::cout << "   v = " << tf << std::endl;
       return tf;
       };
   };

class mVMetricDeriv : public UserFunc<TNT::Vector<double> *, TNT::Vector<double> *>
   {
   protected:
     TNT::Vector<double> cD;
     TubeRegistrator * cRegOp;
   public:
     mVMetricDeriv()
       : cRegOp(0), cD(6, 0.0)
       {
       };
     void use(TubeRegistrator * newRegOp)
       {
       cRegOp = newRegOp;
       };
     TNT::Vector<double> * value(TNT::Vector<double> * m_X)
       {
       cRegOp->SetOffset((*m_X)(1)*offsetUnit, (*m_X)(2)*offsetUnit, (*m_X)(3)*offsetUnit);
       cRegOp->SetAngles((*m_X)(4)*rotUnit, (*m_X)(5)*rotUnit, (*m_X)(6)*rotUnit);
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
       cRegOp->MetricDeriv(&cD(1), &cD(2), &cD(3),
                           &cD(4), &cD(5), &cD(6));
       /*
       std::cout << "  "
                 << cD(1) << " (" << cD(1)/offsetDUnit << ")  "
                 << cD(2) << " (" << cD(2)/offsetDUnit << ")  "
                 << cD(3) << " (" << cD(3)/offsetDUnit << ")  " << std::endl;
       std::cout << "  "
                 << cD(4) << " (" << cD(4)/rotDUnit << ")  "
                 << cD(5) << " (" << cD(5)/rotDUnit << ")  "
                 << cD(6) << " (" << cD(6)/rotDUnit << ")  " << std::endl;
       */
       cD(1) /= offsetDUnit;
       cD(2) /= offsetDUnit;
       cD(3) /= offsetDUnit;
       cD(4) /= rotDUnit;
       cD(5) /= rotDUnit;
       cD(6) /= rotDUnit;
       return &cD;
       };
   };

bool TubeRegistrator::Fit(void)
  {
  mVMetricCost * func = new mVMetricCost;
  func->use(this);

  mVMetricDeriv * funcD = new mVMetricDeriv;
  funcD->use(this);

  OptParabolicFit1D op1D;
  OptimizerND op(6, func, funcD, &op1D);
  op.searchForMin(false);

  double x, y, z, a, b, g;
  x = (*mO)(0);
  y = (*mO)(1);
  z = (*mO)(2);
  a = mAlpha;
  b = mBeta;
  g = mGamma;

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

  cSampling *= 2;
  //cKappa *= 1.5;
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

  cSampling /= 2;
  //cKappa /= 1.5;
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
