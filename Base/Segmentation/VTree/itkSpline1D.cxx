/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkSpline1D.cxx,v $
  Language:  C++
  Date:      $Date: 2003/01/13 19:59:26 $
  Version:   $Revision: 1.3 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkSpline1D.h"
#include <cmath>


namespace itk 
{

class Spline1DValFunc : public UserFunc<double, double> 
{
  Spline1D * spline;
  double cVal;
public:
  Spline1DValFunc(Spline1D * newSpline)
    {
    spline = newSpline;
    cVal = 0;
    };
  const double & value(const double & x)
    {
    cVal = spline->value(x);
    return cVal;
    };
};

class Spline1DDerivFunc : public UserFunc<double, double> 
{
  Spline1D * spline;
  double cDeriv;
public:
  Spline1DDerivFunc(Spline1D * newSpline)
    {
    spline = newSpline;
    cDeriv = 0;
    };
  const double & value(const double & x)
    {
    cDeriv = spline->valueD(x);
    return cDeriv;
    }
};

//
//
//
Spline1D::
Spline1D()
{
  cDefined = false;

  cClip = true;

  cXMin = 0;
  cXMax = 1;

  cOpt1DVal = new Spline1DValFunc(this);
  cOpt1DDeriv = new Spline1DDerivFunc(this);

  use(NULL, NULL);
}

Spline1D::
Spline1D(UserFunc<int, double> *newFuncVal, Optimizer1D *newOpt1D)
{
  cDefined = false;

  cClip = true;

  cXMin = 0;
  cXMax = 1;

  cData.set_size(4);

  cOpt1DVal = new Spline1DValFunc(this);
  cOpt1DDeriv = new Spline1DDerivFunc(this);

  use(newFuncVal, newOpt1D);
}

Spline1D::
~Spline1D()
{
  if(cDefined) 
  {
    cDefined = false;

    delete cOpt1DVal;
    delete cOpt1DDeriv;
  }
}

//
//
//
void Spline1D::
use(UserFunc<int, double> *newFuncVal, Optimizer1D *newOpt1D)
{
    cDefined = true;

    cFuncVal = newFuncVal;

    cOpt1D = newOpt1D;
    if(cOpt1D != NULL)
        cOpt1D->use(cOpt1DVal, cOpt1DDeriv);

    cNewData = true;
}

//
//
//
bool Spline1D::
clipEdge(void)
{
    return cClip;
}

void Spline1D::
clipEdge(bool newClip)
{
    cClip = newClip;
}

int Spline1D::
xMin(void)
{
    return cXMin;
}

void Spline1D::
xMin(int newXMin)
{
    cXMin = newXMin;
}

int Spline1D::
xMax(void)
{
    return cXMax;
}

void Spline1D::
xMax(int newXMax)
{
    cXMax = newXMax;
}

//
//
//
bool Spline1D::
newData(void)
{
    return cNewData;
}

void Spline1D::
newData(bool newNewData)
{
    cNewData = newNewData;
}

//
//
//
void Spline1D::
cGetData(double x)
{
  static int xi = (int)x;

  if((int)x != xi || cNewData) 
    {
    int lx = xi;
    xi = (int)x;
    int offset = xi-lx;
    if( cNewData )
      {
      offset = 100;
      }
    cNewData = false;
    vnl_vector<double> tmpData(4);
    int j=0;
    for(int i=xi-1; i<=xi+2; i++)
      {
      if( j+offset >= 0 && j+offset <= 3 )
        {
        tmpData[j] = cData[j+offset];
        //std::cout << "j = " << j << " : reuse i = " << i 
          //<< " : v = " << tmpData[j] << std::endl;
        j++;
        }
      else
        {
        if(i>=cXMin && i<=cXMax)
          {
          tmpData[j++] = cFuncVal->value(i);
          }
        else
          {
          if(i<cXMin)
            {
            if(cClip)
              {
              tmpData[j++] = 0;
              }
            else
              {
              tmpData[j++] = cFuncVal->value(cXMin) / (1+(cXMin-i));
              }
            }
          else
            {
            if(i>cXMax)
              {
              if(cClip)
                {
                tmpData[j++]= 0;
                }
              else
                {
                tmpData[j++] = cFuncVal->value(cXMax) / (1+(i-cXMax));
                }
              }
            }
          }
        }
      }
    for(j=0; j<4; j++)
      {
      cData[j] = tmpData[j];
      //std::cout << "  cData[" << j << "] = " << cData[j] << std::endl;
      }
    }
}

//
//
//
double Spline1D::
value(double x)
{
    if(!cDefined || (cClip && ((int)x<cXMin || (int)x>cXMax)))
        return 0;

    cGetData(x);

    return dataValue(cData, x - (int)x);
}


double Spline1D::
valueD(double x)
{
    if(!cDefined || (cClip && ((int)x<cXMin || (int)x>cXMax)))
        return 0;

    cGetData(x);

    return dataValueD(cData, x - (int)x);
}


double Spline1D::
valueD2(double x)
{
    if(!cDefined || (cClip && ((int)x<cXMin || (int)x>cXMax)))
        return 0;

    cGetData(x);
        
    return dataValueD2(cData, x - (int)x);
}


double Spline1D::
curv(double x)
{        
    if(!cDefined || (cClip && ((int)x<cXMin || (int)x>cXMax)))
        return 0;

    double xp = valueD(x);
        
    double xpp = valueD2(x);
        
    return xpp/pow(1.0+xp*xp, 1.5);
}

double Spline1D::
valueJet(double x, double * d, double * d2)
{
  if(!cDefined || (cClip && ((int)x<cXMin || (int)x>cXMax)))
  {
    return 0;
  }
  
  cGetData(x);

  return dataValueJet(cData, x - (int)x, d, d2);
}

//
//
//
bool Spline1D::
extreme(double * extX, double * extVal)
{    
    return cOpt1D->extreme(extX, extVal);
}

}; // namespace
