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
#include "tubeSpline1D.h"
#include <cmath>


namespace tube
{

class Spline1DValFunc : public UserFunc<double, double>
{
public:

  Spline1DValFunc(Spline1D * newSpline)
    {
    m_Spline = newSpline;
    m_Val = 0;
    }

  const double & value(const double & x)
    {
    m_Val = m_Spline->value(x);
    return m_Val;
    }

protected:

  Spline1D * m_Spline;
  double     m_Val;

};

class Spline1DDerivFunc : public UserFunc<double, double>
{
public:

  Spline1DDerivFunc(Spline1D * newSpline)
    {
    m_Spline = newSpline;
    m_Deriv = 0;
    }

  const double & value(const double & x)
    {
    m_Deriv = m_Spline->valueD(x);
    return m_Deriv;
    }

protected:

  Spline1D * m_Spline;
  double m_Deriv;

};

//
//
//
Spline1D::
Spline1D()
{
  m_Defined = false;

  m_NewData = true;

  m_Clip = false;

  m_XMin = 0;
  m_XMax = 1;

  m_Data.set_size(4);

  m_Opt1DVal = new Spline1DValFunc(this);
  m_Opt1DDeriv = new Spline1DDerivFunc(this);

  use(NULL, NULL);
}

Spline1D::
Spline1D(UserFunc<int, double> *newFuncVal, Optimizer1D *newOpt1D)
{
  m_Defined = false;

  m_NewData = true;

  m_Clip = false;

  m_XMin = 0;
  m_XMax = 1;

  m_Data.set_size(4);

  m_Opt1DVal = new Spline1DValFunc(this);
  m_Opt1DDeriv = new Spline1DDerivFunc(this);

  use(newFuncVal, newOpt1D);
}

Spline1D::
~Spline1D()
{
  if(m_Defined)
    {
    m_Defined = false;

    delete m_Opt1DVal;
    delete m_Opt1DDeriv;
    }
}

//
//
//
void Spline1D::
use(UserFunc<int, double> *newFuncVal, Optimizer1D *newOpt1D)
{
  m_Defined = true;

  m_FuncVal = newFuncVal;

  m_Opt1D = newOpt1D;
  if(m_Opt1D != NULL)
    {
    m_Opt1D->use(m_Opt1DVal, m_Opt1DDeriv);
    }

  m_NewData = true;
}

//
//
//
bool Spline1D::
clipEdge(void)
{
  return m_Clip;
}

void Spline1D::
clipEdge(bool newClip)
{
  m_Clip = newClip;
}

int Spline1D::
xMin(void)
{
  return m_XMin;
}

void Spline1D::
xMin(int newXMin)
{
  m_XMin = newXMin;
  if(m_Opt1D)
    {
    m_Opt1D->xMin( m_XMin );
    }
}

int Spline1D::
xMax(void)
{
  return m_XMax;
}

void Spline1D::
xMax(int newXMax)
{
  m_XMax = newXMax;
  if(m_Opt1D)
    {
    m_Opt1D->xMax( m_XMax );
    }
}

//
//
//
bool Spline1D::
newData(void)
{
  return m_NewData;
}

void Spline1D::
newData(bool newNewData)
{
  m_NewData = newNewData;
}

//
//
//
void Spline1D::
m_GetData(double x)
{
  static int xi = (int)x;

  if((int)x != xi || m_NewData)
    {
    int lx = xi;
    xi = (int)x;
    int offset = xi-lx;
    if( m_NewData )
      {
      offset = 100;
      }
    m_NewData = false;
    vnl_vector<double> tmpData(4);
    int j=0;
    for(int i=xi-1; i<=xi+2; i++)
      {
      if( j+offset >= 0 && j+offset <= 3 )
        {
        tmpData[j] = m_Data[j+offset];
        j++;
        }
      else
        {
        if(i>=m_XMin && i<=m_XMax)
          {
          tmpData[j++] = m_FuncVal->value(i);
          }
        else
          {
          if(i<m_XMin)
            {
            if(m_Clip)
              {
              tmpData[j++] = m_FuncVal->value(m_XMin);
              }
            else
              {
              int tmpI = m_XMin + (m_XMin-i);
              if( tmpI > m_XMax )
                {
                tmpI = m_XMax;
                }
              tmpData[j++] = m_FuncVal->value(tmpI);
              }
            }
          else
            {
            if(i>m_XMax)
              {
              if(m_Clip)
                {
                tmpData[j++] = m_FuncVal->value(m_XMax);
                }
              else
                {
                int tmpI = m_XMax - (i-m_XMax);
                if( tmpI < m_XMin )
                  {
                  tmpI = m_XMin;
                  }
                tmpData[j++] = m_FuncVal->value(tmpI);
                }
              }
            }
          }
        }
      }
    for(j=0; j<4; j++)
      {
      m_Data[j] = tmpData[j];
      //std::cout << "  m_Data[" << j << "] = " << m_Data[j] << std::endl;
      }
    }
}

//
//
//
double Spline1D::
value(double x)
{
  if(!m_Defined || (m_Clip && ((int)x<m_XMin || (int)x>m_XMax)))
  //if(!m_Defined)
    {
    return 0;
    }

  m_GetData(x);

  return dataValue(m_Data, x - (int)x);
}


double Spline1D::
valueD(double x)
{
  if(!m_Defined || (m_Clip && ((int)x<m_XMin || (int)x>m_XMax)))
  //if(!m_Defined)
    {
    return 0;
    }

  m_GetData(x);

  return dataValueD(m_Data, x - (int)x);
}


double Spline1D::
valueD2(double x)
{
  if(!m_Defined || (m_Clip && ((int)x<m_XMin || (int)x>m_XMax)))
  //if(!m_Defined)
    {
    return 0;
    }

  m_GetData(x);

  return dataValueD2(m_Data, x - (int)x);
}


double Spline1D::
curv(double x)
{
  if(!m_Defined || (m_Clip && ((int)x<m_XMin || (int)x>m_XMax)))
  //if(!m_Defined)
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
  if(!m_Defined || (m_Clip && ((int)x<m_XMin || (int)x>m_XMax)))
  //if(!m_Defined)
    {
    return 0;
    }

  m_GetData(x);

  return dataValueJet(m_Data, x - (int)x, d, d2);
}

//
//
//
bool Spline1D::
extreme(double * extX, double * extVal)
{
  return m_Opt1D->extreme(extX, extVal);
}

//
//
//
void Spline1D::
PrintSelf( std::ostream & os ) const
{
  if( m_Defined )
    {
    os << "m_Defined = True" << std::endl;
    }
  else
    {
    os << "m_Defined = False" << std::endl;
    }
  os << "m_FuncVal = " << m_FuncVal << std::endl;
  if( m_Clip )
    {
    os << "m_Clip = True" << std::endl;
    }
  else
    {
    os << "m_Clip = False" << std::endl;
    }
  os << "m_XMin = " << m_XMin << std::endl;
  os << "m_XMax = " << m_XMax << std::endl;
  if( m_NewData )
    {
    os << "m_NewData = True" << std::endl;
    }
  else
    {
    os << "m_NewData = False" << std::endl;
    }
  os << "m_Data = " << m_Data << std::endl;
  os << "m_Opt1DVal = " << m_Opt1DVal << std::endl;
  os << "m_Opt1DDeriv = " << m_Opt1DDeriv << std::endl;
  if( m_Opt1D != NULL )
    {
    os << "m_Opt1D = " << std::endl;
    m_Opt1D->PrintSelf( os );
    }
  else
    {
    os << "m_Opt1D = NULL" << std::endl;
    }
}

} // namespace
