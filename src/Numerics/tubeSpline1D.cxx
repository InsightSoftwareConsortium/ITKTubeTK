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

#include "tubeSpline1D.h"

namespace tube
{

class Spline1DValueFunction : public UserFunction< double, double >
{
public:

  typedef Spline1DValueFunction           Self;
  typedef UserFunction< double, double >  Superclass;
  typedef Self *                          Pointer;
  typedef const Self *                    ConstPointer;

  typedef double                          InputType;
  typedef double                          OutputType;

  tubeTypeMacro( Spline1DValueFunction );

  Spline1DValueFunction( Spline1D::Pointer spline )
    {
    m_Spline = spline;
    m_Value = 0.0;
    }

  const OutputType & Value( const InputType & input ) override
    {
    m_Value = m_Spline->Value( input );
    return m_Value;
    }

private:

  Spline1DValueFunction( const Self & self );
  void operator=( const Self & self );

  Spline1D::Pointer  m_Spline;
  OutputType         m_Value;

}; // End class Spline1DValueFunction

class Spline1DDerivativeFunction : public UserFunction< double, double >
{
public:

  typedef Spline1DDerivativeFunction      Self;
  typedef UserFunction< double, double >  Superclass;
  typedef Self *                          Pointer;
  typedef const Self *                    ConstPointer;

  typedef double                          InputType;
  typedef double                          OutputType;

  tubeTypeMacro( Spline1DDerivativeFunction );

  Spline1DDerivativeFunction( Spline1D::Pointer spline )
    {
    m_Spline = spline;
    m_Derivative = 0.0;
    }

  const OutputType & Value( const InputType & input ) override
    {
    m_Derivative = m_Spline->ValueD( input );
    return m_Derivative;
    }

private:

  Spline1DDerivativeFunction( const Self & self );
  void operator=( const Self & self );

  Spline1D::Pointer  m_Spline;
  OutputType         m_Derivative;

}; // End class Spline1DDerivativeFunction


Spline1D
::Spline1D( void )
  : m_Data( 4, 0.0 )
{
  m_Defined = false;

  m_NewData = true;

  m_Clip = false;

  m_XMin = 0;
  m_XMax = 1;

  m_Optimizer1DVal = new Spline1DValueFunction( this );
  m_Optimizer1DDeriv = new Spline1DDerivativeFunction( this );

  m_Optimizer1D = NULL;
  m_FuncVal = NULL;

  this->Use( m_FuncVal, m_Optimizer1D );
}


Spline1D
::Spline1D( UserFunction< int, double >::Pointer funcVal,
  Optimizer1D::Pointer optimizer1D )
  : m_Data( 4, 0.0 )
{
  m_Defined = false;

  m_NewData = true;

  m_Clip = false;

  m_XMin = 0;
  m_XMax = 1;

  m_Optimizer1DVal = new Spline1DValueFunction( this );
  m_Optimizer1DDeriv = new Spline1DDerivativeFunction( this );

  m_Optimizer1D = optimizer1D;
  m_FuncVal = funcVal;

  this->Use( m_FuncVal, m_Optimizer1D );
}


Spline1D
::~Spline1D( void )
{
  if( m_Defined )
    {
    m_Defined = false;

    delete m_Optimizer1DVal;
    delete m_Optimizer1DDeriv;
    }
}


void
Spline1D
::Use( ValueFunctionType::Pointer funcVal, Optimizer1D::Pointer optimizer1D )
{
  m_Defined = true;

  m_FuncVal = funcVal;

  m_Optimizer1D = optimizer1D;
  if( m_Optimizer1D != NULL )
    {
    m_Optimizer1D->Use( m_Optimizer1DVal, m_Optimizer1DDeriv );
    }

  m_NewData = true;
}


void
Spline1D
::SetXMin( int xMin )
{
  m_XMin = xMin;

  if( m_Optimizer1D )
    {
    m_Optimizer1D->SetXMin( m_XMin );
    }
}


void
Spline1D
::SetXMax( int xMax )
{
  m_XMax = xMax;

  if( m_Optimizer1D )
    {
    m_Optimizer1D->SetXMax( m_XMax );
    }
}


void
Spline1D
::m_GetData( double x )
{
  static int xi = ( int )x;

  if( ( int )x != xi || m_NewData )
    {
    int lx = xi;
    xi = ( int )x;
    int offset = xi-lx;
    if( m_NewData )
      {
      offset = 100;
      }
    m_NewData = false;
    vnl_vector<double> tmpData( 4 );
    int j=0;
    for( int i=xi-1; i<=xi+2; i++ )
      {
      if( j+offset >= 0 && j+offset <= 3 )
        {
        tmpData[j] = m_Data[j+offset];
        j++;
        }
      else
        {
        if( i>=m_XMin && i<=m_XMax )
          {
          tmpData[j++] = m_FuncVal->Value( i );
          }
        else
          {
          if( i<m_XMin )
            {
            if( m_Clip )
              {
              tmpData[j++] = m_FuncVal->Value( m_XMin );
              }
            else
              {
              int tmpI = m_XMin + ( m_XMin-i );
              if( tmpI > m_XMax )
                {
                tmpI = m_XMax;
                }
              tmpData[j++] = m_FuncVal->Value( tmpI );
              }
            }
          else
            {
            if( i>m_XMax )
              {
              if( m_Clip )
                {
                tmpData[j++] = m_FuncVal->Value( m_XMax );
                }
              else
                {
                int tmpI = m_XMax - ( i-m_XMax );
                if( tmpI < m_XMin )
                  {
                  tmpI = m_XMin;
                  }
                tmpData[j++] = m_FuncVal->Value( tmpI );
                }
              }
            }
          }
        }
      }
    for( j=0; j<4; j++ )
      {
      m_Data[j] = tmpData[j];
      }
    }
}


double
Spline1D
::Value( double x )
{
  if( !m_Defined || ( m_Clip && ( x<( double )m_XMin || x>( double )m_XMax ) ) )
    {
    std::cout << "clipping: " << m_XMin << " <= " << x << " <= " << m_XMax
      << std::endl;
    return 0;
    }

  this->m_GetData( x );

  return this->DataValue( m_Data, x - ( int )x );
}


double
Spline1D
::ValueD( double x )
{
  if( !m_Defined || ( m_Clip && ( x<( double )m_XMin || x>( double )m_XMax ) ) )
    {
    return 0;
    }

  this->m_GetData( x );

  return this->DataValueD( m_Data, x - ( int )x );
}


double
Spline1D
::ValueD2( double x )
{
  if( !m_Defined || ( m_Clip && ( x<( double )m_XMin || x>( double )m_XMax ) ) )
    {
    return 0;
    }

  this->m_GetData( x );

  return this->DataValueD2( m_Data, x - ( int )x );
}


double
Spline1D
::Curv( double x )
{
  if( !m_Defined || ( m_Clip && ( x<( double )m_XMin || x>( double )m_XMax ) ) )
    {
    return 0;
    }

  double xp = ValueD( x );

  double xpp = ValueD2( x );

  return xpp/std::pow( 1.0+xp*xp, 1.5 );
}


double
Spline1D
::ValueJet( double x, double * d, double * d2 )
{
  if( !m_Defined || ( m_Clip && ( x<( double )m_XMin || x>( double )m_XMax ) ) )
    {
    return 0;
    }

  this->m_GetData( x );

  return this->DataValueJet( m_Data, x - ( int )x, d, d2 );
}


bool
Spline1D
::Extreme( double * extX, double * extVal )
{
  return m_Optimizer1D->Extreme( extX, extVal );
}


void
Spline1D
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  os << indent << "Defined:          " << m_Defined << std::endl;
  os << indent << "FuncVal:          " << m_FuncVal << std::endl;
  os << indent << "Clip:             " << m_Clip << std::endl;
  os << indent << "XMin:             " << m_XMin << std::endl;
  os << indent << "XMax:             " << m_XMax << std::endl;
  os << indent << "NewData:          " << m_NewData << std::endl;
  os << indent << "Data:             " << m_Data << std::endl;
  os << indent << "Optimizer1DVal:   " << m_Optimizer1DVal << std::endl;
  os << indent << "Optimizer1DDeriv: " << m_Optimizer1DDeriv << std::endl;
  os << indent << "Optimizer1D:      " << m_Optimizer1D << std::endl;
}

} // End namespace tube
