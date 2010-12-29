/*=========================================================================

Library:   TubeTK/VTree

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
#include "tubeSpline1D.h"
#include <cmath>


namespace tube
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

  cClip = false;

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

  cClip = false;

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
      {
      cOpt1D->use(cOpt1DVal, cOpt1DDeriv);
      }

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
    if(cOpt1D)
      {
      cOpt1D->xMin( cXMin );
      }
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
    if(cOpt1D)
      {
      cOpt1D->xMax( cXMax );
      }
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
              tmpData[j++] = cFuncVal->value(cXMin);
              }
            else
              {
              int tmpI = cXMin + (cXMin-i);
              if( tmpI > cXMax )
                {
                tmpI = cXMax;
                }
              tmpData[j++] = cFuncVal->value(tmpI);
              }
            }
          else
            {
            if(i>cXMax)
              {
              if(cClip)
                {
                tmpData[j++] = cFuncVal->value(cXMax);
                }
              else
                {
                int tmpI = cXMax - (i-cXMax);
                if( tmpI < cXMin )
                  {
                  tmpI = cXMin;
                  }
                tmpData[j++] = cFuncVal->value(tmpI);
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
  //if(!cDefined)
    {
    return 0;
    }

  cGetData(x);

  return dataValue(cData, x - (int)x);
}


double Spline1D::
valueD(double x)
{
  if(!cDefined || (cClip && ((int)x<cXMin || (int)x>cXMax)))
  //if(!cDefined)
    {
    return 0;
    }

  cGetData(x);

  return dataValueD(cData, x - (int)x);
}


double Spline1D::
valueD2(double x)
{
  if(!cDefined || (cClip && ((int)x<cXMin || (int)x>cXMax)))
  //if(!cDefined)
    {
    return 0;
    }

  cGetData(x);

  return dataValueD2(cData, x - (int)x);
}


double Spline1D::
curv(double x)
{
  if(!cDefined || (cClip && ((int)x<cXMin || (int)x>cXMax)))
  //if(!cDefined)
    {
    return 0;
    }

  double xp = valueD(x);

  double xpp = valueD2(x);

  return xpp/pow(1.0+xp*xp, 1.5);
}

double Spline1D::
valueJet(double x, double * d, double * d2)
{
  if(!cDefined || (cClip && ((int)x<cXMin || (int)x>cXMax)))
  //if(!cDefined)
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

//
//
//
void Spline1D::
PrintSelf( std::ostream & os ) const
{
  if( cDefined )
    {
    os << "cDefined = True" << std::endl;
    }
  else
    {
    os << "cDefined = False" << std::endl;
    }
  os << "cFuncVal = " << cFuncVal << std::endl;
  if( cClip )
    {
    os << "cClip = True" << std::endl;
    }
  else
    {
    os << "cClip = False" << std::endl;
    }
  os << "cXMin = " << cXMin << std::endl;
  os << "cXMax = " << cXMax << std::endl;
  if( cNewData )
    {
    os << "cNewData = True" << std::endl;
    }
  else
    {
    os << "cNewData = False" << std::endl;
    }
  os << "cData = " << cData << std::endl;
  os << "cOpt1DVal = " << cOpt1DVal << std::endl;
  os << "cOpt1DDeriv = " << cOpt1DDeriv << std::endl;
  if( cOpt1D != NULL )
    {
    os << "cOpt1D = " << std::endl;
    cOpt1D->PrintSelf( os );
    }
  else
    {
    os << "cOpt1D = NULL" << std::endl;
    }
}

}; // namespace
