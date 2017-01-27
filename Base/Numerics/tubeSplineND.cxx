/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include "tubeSplineND.h"

#include <itkImageRegionIterator.h>

namespace tube
{

class SplineNDValueFunction
: public UserFunction< vnl_vector< double >, double >
{
public:

  typedef SplineNDValueFunction                         Self;
  typedef UserFunction< vnl_vector< double >, double >  Superclass;
  typedef Self *                                        Pointer;
  typedef const Self *                                  ConstPointer;

  typedef SplineND::VectorType                          VectorType;
  typedef VectorType                                    InputType;
  typedef double                                        OutputType;

  tubeTypeMacro( SplineNDValueFunction );

  SplineNDValueFunction( SplineND::Pointer spline )
    {
    m_Spline = spline;
    m_Value = 0.0;
    }

  const OutputType & Value( const InputType & input )
    {
    m_Value = m_Spline->Value( input );
    return m_Value;
    }

private:

  SplineNDValueFunction( const Self & self );
  void operator=( const Self & self );

  SplineND::Pointer  m_Spline;
  OutputType         m_Value;

}; // End class SplineNDValueFunction

class SplineNDDerivativeFunction
  : public UserFunction< vnl_vector< double >, vnl_vector< double > >
{
public:

  typedef SplineNDDerivativeFunction  Self;
  typedef UserFunction< vnl_vector< double >, vnl_vector< double > >
                                      Superclass;
  typedef Self *                      Pointer;
  typedef const Self *                ConstPointer;

  typedef SplineND::VectorType        VectorType;
  typedef VectorType                  InputType;
  typedef VectorType                  OutputType;

  tubeTypeMacro( SplineNDDerivativeFunction );

  SplineNDDerivativeFunction( SplineND::Pointer spline )
    {
    m_Spline = spline;
    m_Derivative.set_size( m_Spline->GetDimension() );
    }

  const OutputType & Value( const InputType & input )
    {
    m_Derivative = m_Spline->ValueD( input );
    return m_Derivative;
    }

private:

  SplineNDDerivativeFunction( const Self & self );
  void operator=( const Self & self );

  SplineND::Pointer  m_Spline;
  OutputType         m_Derivative;

}; // End class SplineNDDerivativeFunction


SplineND
::SplineND( void )
  : m_FuncVal( NULL )
{
  m_Dimension = 0;

  m_Clip = false;
  m_NewData = true;

  m_Val = 0;
  m_Data = NULL;
  m_DataWS = NULL;
  m_DataWSX = NULL;
  m_DataWSXX = NULL;

  m_OptimizerNDVal = new SplineNDValueFunction( this );
  m_OptimizerNDDeriv = new SplineNDDerivativeFunction( this );

  m_OptimizerND = NULL;
  m_Spline1D = NULL;

  this->Use( 0, NULL, NULL, NULL );
}


SplineND
::SplineND( unsigned int dimension,
  ValueFunctionType::Pointer funcVal,
  Spline1D::Pointer spline1D,
  Optimizer1D::Pointer optimizer1D )
{
  m_Dimension = 0;

  m_Clip = false;
  m_NewData = true;

  m_Val = 0;
  m_Data = NULL;
  m_DataWS = NULL;
  m_DataWSX = NULL;
  m_DataWSXX = NULL;

  m_FuncVal = NULL;

  m_OptimizerNDVal = new SplineNDValueFunction( this );
  m_OptimizerNDDeriv = new SplineNDDerivativeFunction( this );

  m_OptimizerND = NULL;
  m_Spline1D = NULL;

  this->Use( dimension, funcVal, spline1D, optimizer1D );
}


SplineND
::~SplineND( void )
{
  delete m_OptimizerNDVal;
  delete m_OptimizerNDDeriv;
  if( m_OptimizerND != NULL )
    {
    delete m_OptimizerND;
    }
}


void
SplineND
::Use( unsigned int dimension,
  ValueFunctionType::Pointer funcVal,
  Spline1D::Pointer spline1D,
  Optimizer1D::Pointer optimizer1D )
{
  if( m_OptimizerND != NULL )
    {
    delete m_OptimizerND;
    m_OptimizerND = NULL;
    }

  m_Dimension = dimension;
  m_XMin.set_size( m_Dimension );
  m_XMin.fill( ( int )0 );
  m_XMax.set_size( m_Dimension );
  m_XMax.fill( ( int )1 );
  m_Xi.set_size( m_Dimension );
  m_D.set_size( m_Dimension );
  m_H.set_size( m_Dimension, m_Dimension );

  ImageType::SizeType dimSize;
  for( unsigned int i=0; i<dimSize.GetSizeDimension(); i++ )
    {
    dimSize[i] = 1;
    }
  for( unsigned int i=0; i<m_Dimension; i++ )
    {
    dimSize[i] = 4;
    }

  m_Data   = ImageType::New();
  ImageType::RegionType region;
  region.SetSize( dimSize );
  m_Data->SetRegions( region );
  m_Data->Allocate();

  m_DataWS = ImageType::New();
  m_DataWS->SetRegions( region );
  m_DataWS->Allocate();

  m_DataWSX = VectorImageType::New();
  m_DataWSXX = VectorImageType::New();
  m_DataWSX->Reserve( m_Dimension );
  m_DataWSXX->Reserve( m_Dimension );

  VectorImageType::Iterator itWSX  = m_DataWSX->Begin();
  VectorImageType::Iterator itWSXX = m_DataWSXX->Begin();
  while( itWSX != m_DataWSX->End() )
    {
    itWSX->Value() = ImageType::New();

    itWSX->Value()->SetLargestPossibleRegion( region );
    itWSX->Value()->SetBufferedRegion( region );
    itWSX->Value()->SetRequestedRegion( region );
    itWSX->Value()->Allocate();

    itWSXX->Value() = ImageType::New();
    itWSXX->Value()->SetLargestPossibleRegion( region );
    itWSXX->Value()->SetBufferedRegion( region );
    itWSXX->Value()->SetRequestedRegion( region );
    itWSXX->Value()->Allocate();

    itWSX++;
    itWSXX++;
    }

  m_Data1D.set_size( 4 );

  m_Val = 0;
  m_FuncVal = funcVal;
  m_Spline1D = spline1D;

  if( optimizer1D != NULL )
    {
    m_OptimizerND = new OptimizerND( m_Dimension,
      m_OptimizerNDVal,
      m_OptimizerNDDeriv,
      optimizer1D );
    }

  m_NewData = true;
}


void
SplineND
::SetXMin( IntVectorType xMin )
{
  VectorType t( m_Dimension );

  for( unsigned int i = 0; i < m_Dimension; ++i )
    {
    t( i ) = xMin( i );
    }

  m_OptimizerND->SetXMin( t );
  m_XMin = xMin;
}


void
SplineND
::SetXMax( IntVectorType xMax )
{
  VectorType t( m_Dimension );

  for( unsigned int i = 0; i < m_Dimension; ++i )
    {
    t( i ) = xMax( i );
    }

  m_OptimizerND->SetXMax( t );
  m_XMax = xMax;
}


void
SplineND
::m_GetData( const VectorType & x )
{
  bool eql = true;
  if( !m_NewData )
    {
    for( unsigned int i=0; i<m_Dimension; i++ )
      {
      if( m_Xi( i ) != ( int )x( i ) )
        {
        eql = false;
        break;
        }
      }
    }
  if( !eql || m_NewData )
    {
    if( m_NewData )
      {
      m_NewData = false;
      for( unsigned int i=0; i<m_Dimension; i++ )
        {
        m_Xi( i ) = ( int )x( i );
        }

      itk::ImageRegionIterator<ImageType> it( m_Data,
        m_Data->GetLargestPossibleRegion() );
      it.GoToBegin();

      IntVectorType p( m_Dimension );
      IntVectorType xiOffset( m_Dimension, -1 );
      bool done = false;
      while( !done )
        {
        p = m_Xi + xiOffset;
        for( unsigned int i=0; i<m_Dimension; i++ )
          {
          if( p( i ) < m_XMin( i ) )
            {
            if( m_Clip )
              {
              p( i ) = m_XMin( i );
              }
            else
              {
              p( i ) = m_XMin( i ) + ( m_XMin( i ) - p( i ) );
              if( p( i ) > m_XMax( i ) )
                {
                p( i ) = m_XMax( i );
                }
              }
            }
          else if( p( i ) > m_XMax( i ) )
            {
            if( m_Clip )
              {
              p( i ) = m_XMax( i );
              }
            else
              {
              p( i ) = m_XMax( i ) - ( p( i ) - m_XMax( i ) );
              if( p( i ) < m_XMin( i ) )
                {
                p( i ) = m_XMin( i );
                }
              }
            }
          }
        it.Set( m_FuncVal->Value( p ) );
        ++it;
        unsigned int dim = 0;
        while( !done && dim<m_Dimension && ( ++xiOffset( dim ) )>2 )
          {
          xiOffset( dim++ ) = -1;
          if( dim >= m_Dimension )
            {
            done = true;
            }
          }
        }
      }
    else
      {
      IntVectorType p( m_Dimension );
      IntVectorType pOld( m_Dimension );
      IntVectorType xiOffset( m_Dimension );
      IntVectorType shift( m_Dimension, 100.0 );

      for( unsigned int i=0; i<m_Dimension; i++ )
        {
        shift[i] = ( int )x( i ) - m_Xi( i );
        }

      for( unsigned int i=0; i<m_Dimension; i++ )
        {
        m_Xi( i ) = ( int )x( i );
        }

      itk::ImageRegionIterator<ImageType> it( m_DataWS,
        m_DataWS->GetLargestPossibleRegion() );
      it.GoToBegin();

      xiOffset = -1;
      bool done = false;
      while( !done )
        {
        p = m_Xi + xiOffset;
        pOld = xiOffset + shift + 1;
        for( unsigned int i=0; i<m_Dimension; i++ )
          {
          if( p( i ) < m_XMin( i ) )
            {
            if( m_Clip )
              {
              p( i ) = m_XMin( i );
              }
            else
              {
              p( i ) = m_XMin( i ) + ( m_XMin( i ) - p( i ) );
              if( p( i ) > m_XMax( i ) )
                {
                p( i ) = m_XMax( i );
                }
              }
            }
          else if( p( i ) > m_XMax( i ) )
            {
            if( m_Clip )
              {
              p( i ) = m_XMax( i );
              }
            else
              {
              p( i ) = m_XMax( i ) - ( p( i ) - m_XMax( i ) );
              if( p( i ) < m_XMin( i ) )
                {
                p( i ) = m_XMin( i );
                }
              }
            }
          }
        bool reuse = true;
        for( unsigned int j=0; j<m_Dimension; j++ )
          {
          if( pOld[j] < 0 || pOld[j] > 3 )
            {
            reuse = false;
            break;
            }
          }
        if( reuse )
          {
          ImageType::IndexType indx;
          indx.Fill( 0 );
          for( unsigned int i=0; i<m_Dimension; i++ )
            {
            indx[i] = pOld[i];
            }
          it.Set( m_Data->GetPixel( indx ) );
          }
        else
          {
          it.Set( m_FuncVal->Value( p ) );
          }
        ++it;
        unsigned int dim = 0;
        while( !done && dim<m_Dimension && ( ++xiOffset( dim ) )>2 )
          {
          xiOffset( dim++ ) = -1;
          if( dim >= m_Dimension )
            {
            done = true;
            }
          }
        }

      itk::ImageRegionIterator<ImageType> itDest( m_Data,
        m_Data->GetLargestPossibleRegion() );
      it.GoToBegin();
      itDest.GoToBegin();
      while( !it.IsAtEnd() )
        {
        itDest.Set( it.Get() );
        ++it;
        ++itDest;
        }
      }
    }
}


OptimizerND ::Pointer
SplineND
::GetOptimizerND( void )
{
  return m_OptimizerND;
}


double
SplineND
::Value( const VectorType & x )
{
  this->m_GetData( x );

  itk::ImageRegionIterator<ImageType> itData( m_Data,
    m_Data->GetLargestPossibleRegion() );

  itk::ImageRegionIterator<ImageType> itDataWS( m_DataWS,
    m_DataWS->GetLargestPossibleRegion() );

  itData.GoToBegin();
  itDataWS.GoToBegin();
  while( !itData.IsAtEnd() )
    {
    itDataWS.Set( itData.Get() );
    ++itDataWS;
    ++itData;
    }

  itk::ImageRegionIterator<ImageType> itDataWSColumn( m_DataWS,
    m_DataWS->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> itDataWSDest( m_DataWS,
    m_DataWS->GetLargestPossibleRegion() );

  for( int i=( int )m_Dimension-1; i>=0; i-- )
    {
    unsigned int k = ( unsigned int )std::pow( ( float )4, ( int )i );
    itDataWSColumn.GoToBegin();
    itDataWSDest.GoToBegin();
    for( unsigned int j=0; j<k; j++ )
      {
      for( unsigned int ind=0; ind<4; ind++ )
        {
        m_Data1D( ind )= itDataWSColumn.Get();
        ++itDataWSColumn;
        }

      itDataWSDest.Set( m_Spline1D->DataValue( m_Data1D,
        ( ( x( ( int )m_Dimension-i-1 )
        - ( int )x( ( int )m_Dimension-i-1 ) ) ) ) );
      ++itDataWSDest;
      }
    }

  itDataWS.GoToBegin();
  m_Val = itDataWS.Get();
  return m_Val;
}


double
SplineND
::ValueD( const VectorType & x, IntVectorType & dx )
{
  this->m_GetData( x );

  itk::ImageRegionIterator<ImageType> itData( m_Data,
    m_Data->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> itDataWS( m_DataWS,
    m_DataWS->GetLargestPossibleRegion() );

  itData.GoToBegin();
  itDataWS.GoToBegin();
  while( !itData.IsAtEnd() )
    {
    itDataWS.Set( itData.Get() );
    ++itDataWS;
    ++itData;
    }

  itk::ImageRegionIterator<ImageType> itDataWSColumn( m_DataWS,
    m_DataWS->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> itDataWSDest( m_DataWS,
    m_DataWS->GetLargestPossibleRegion() );

  for( int i=( int )m_Dimension-1; i>=0; i-- )
    {
    itDataWSColumn.GoToBegin();
    itDataWSDest.GoToBegin();
    unsigned int k = ( unsigned int )std::pow( ( float )4, ( int )i );
    switch( dx( ( int )m_Dimension-i-1 ) )
      {
      default:
      case 0:
        for( unsigned int j=0; j<k; j++ )
          {
          for( unsigned int ind=0; ind<4; ind++ )
            {
            m_Data1D( ind )= itDataWSColumn.Get();
            ++itDataWSColumn;
            }

          itDataWSDest.Set( m_Spline1D->DataValue( m_Data1D,
            ( ( x( ( int )m_Dimension-i-1 ) 
            - ( int )x( ( int )m_Dimension-i-1 ) ) ) ) );
          ++itDataWSDest;
          }
        break;
      case 1:
        for( unsigned int j=0; j<k; j++ )
          {
          for( unsigned int ind=0; ind<4; ind++ )
            {
            m_Data1D( ind )= itDataWSColumn.Get();
            ++itDataWSColumn;
            }

          itDataWSDest.Set( m_Spline1D->DataValueD( m_Data1D,
            ( ( x( ( int )m_Dimension-i-1 )
            - ( int )x( ( int )m_Dimension-i-1 ) ) ) ) );
          ++itDataWSDest;
          }
        break;
      case 2:
        for( unsigned int j=0; j<k; j++ )
          {
          for( unsigned int ind=0; ind<4; ind++ )
            {
            m_Data1D( ind )= itDataWSColumn.Get();
            ++itDataWSColumn;
            }

          itDataWSDest.Set( m_Spline1D->DataValueD2( m_Data1D,
            ( ( x( ( int )m_Dimension-i-1 ) 
            - ( int )x( ( int )m_Dimension-i-1 ) ) ) ) );
          ++itDataWSDest;
          }
        break;
      }
    }

  itDataWS.GoToBegin();
  m_Val = itDataWS.Get();
  return m_Val;
}

SplineND::VectorType &
SplineND
::ValueD( const VectorType & x )
{
  IntVectorType dx( m_Dimension, 0 );

  for( unsigned int i=0; i<m_Dimension; i++ )
    {
    dx( i ) = 1;
    m_D( i ) = ValueD( x, dx );
    dx( i ) = 0;
    }

  return m_D;
}


SplineND::MatrixType &
SplineND
::Hessian( const VectorType & x )
{

  IntVectorType dx( m_Dimension, 0 );

  VectorType d( m_Dimension );
  VectorType d2( m_Dimension );
  ValueVDD2( x, d, d2 );
  for( unsigned int i=0; i<m_Dimension; i++ )
    {
    m_H.put( i, i, d2( i ) );
    m_D( i ) = d( i );
    }

  for( unsigned int i=0; i<m_Dimension; i++ )
    {
    for( unsigned int j=i+1; j<m_Dimension; j++ )
      {
      dx( i ) = 1;
      dx( j ) = 1;
      m_H.put( i, j, ValueD( x, dx ) );
      m_H.put( j, i, m_H.get( i, j ) );
      dx( i ) = 0;
      dx( j ) = 0;
      }
    }
  return m_H;
}


double
SplineND
::ValueJet( const VectorType & x, VectorType & d, MatrixType & h )
{
  IntVectorType dx( m_Dimension, 0 );

  VectorType lD( m_Dimension );
  VectorType lD2( m_Dimension );
  double v = ValueVDD2( x, lD, lD2 );
  MatrixType tempMatrix( m_Dimension, m_Dimension );

  for( unsigned int i=0; i< m_Dimension; i++ )
    {
    m_H.put( i, i, lD2( i ) );
    m_D( i ) = lD( i );
    }

  for( unsigned int i=0; i<m_Dimension; i++ )
    {
    for( unsigned int j=i+1; j<m_Dimension; j++ )
      {
      dx( i ) = 1;
      dx( j ) = 1;
      m_H.put( i, j, ValueD( x, dx ) );
      m_H.put( j, i, m_H.get( i, j ) );
      dx( i ) = 0;
      dx( j ) = 0;
      }
    }

  for( unsigned int i=0; i<m_Dimension; i++ )
    {
    d( i ) = m_D( i );
    }

  h = m_H;

  return v;
}


double
SplineND
::ValueVDD2( const VectorType & x, VectorType & d, VectorType & d2 )
{
  this->m_GetData( x );

  VectorImageType::Iterator itWSX  = m_DataWSX->Begin();
  VectorImageType::Iterator itWSXX = m_DataWSXX->Begin();

  itk::ImageRegionIterator<ImageType> itData( m_Data,
    m_Data->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> itDataWS( m_DataWS,
    m_DataWS->GetLargestPossibleRegion() );

  /** Could optimize that */
  itDataWS.GoToBegin();
  itData.GoToBegin();
  while( !itData.IsAtEnd() )
    {
    itDataWS.Set( itData.Get() );
    ++itDataWS;
    ++itData;
    }

  while( itWSX != m_DataWSX->End() )
    {
    /** Iterators to copy image information */
    itk::ImageRegionIterator<ImageType> itImageWSX( itWSX->Value(),
      itWSX->Value()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<ImageType> itImageWSXX( itWSXX->Value(),
      itWSXX->Value()->GetLargestPossibleRegion() );
    itData.GoToBegin();
    itImageWSX.GoToBegin();
    itImageWSXX.GoToBegin();
    while( !itImageWSX.IsAtEnd() )
      {
      itImageWSX.Set( itData.Get() );
      itImageWSXX.Set( itData.Get() );
      ++itImageWSX;
      ++itImageWSXX;
      ++itData;
      }
    itWSX++;
    itWSXX++;
    }

  itWSX  = m_DataWSX->Begin();
  itWSXX = m_DataWSXX->Begin();
  double vD;
  double vD2;

  for( int i=( int )m_Dimension-1; i>=0; i-- )
    {
    unsigned int k = ( unsigned int )std::pow( ( float )4, ( int )i );

    itk::ImageRegionIterator<ImageType> itImageWSX( itWSX->Value(),
      itWSX->Value()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<ImageType> itImageWSXX( itWSXX->Value(),
      itWSXX->Value()->GetLargestPossibleRegion() );
    itImageWSX.GoToBegin();
    itImageWSXX.GoToBegin();

    itk::ImageRegionIterator<ImageType> itDataWSX( m_DataWS,
      m_DataWS->GetLargestPossibleRegion() );
    itDataWSX.GoToBegin();

    for( unsigned int j=0; j<k; j++ )
      {
      itDataWSX.GoToBegin();
      for( unsigned int offset=0; offset<j*4; offset++ )
        {
        ++itDataWSX;
        }
      for( unsigned int ind=0; ind<4; ind++ )
        {
        m_Data1D( ind )= itDataWSX.Get();
        ++itDataWSX;
        }

      itDataWSX.GoToBegin();
      for( unsigned int offset=0; offset<j; offset++ )
        {
        ++itDataWSX;
        }

      itDataWSX.Set( m_Spline1D->DataValueJet( m_Data1D,
        ( ( x( ( int )m_Dimension-i-1 ) 
        - ( int )x( ( int )m_Dimension-i-1 ) ) ), &vD, &vD2 ) );

      for( unsigned int l=0; l<m_Dimension; l++ )
        {

        if( ( int )m_Dimension-i != ( int )l )
          {
          VectorImageType::Iterator itWSX2  = m_DataWSX->Begin();
          VectorImageType::Iterator itWSXX2 = m_DataWSXX->Begin();

          for( unsigned int offset=0; offset<l; offset++ )
            {
            itWSX2++;
            itWSXX2++;
            }
          itk::ImageRegionIterator<ImageType> itImageWSX2( itWSX2->Value(),
            itWSX2->Value()->GetLargestPossibleRegion() );
          itk::ImageRegionIterator<ImageType> itImageWSXX2(
            itWSXX2->Value(),
            itWSXX2->Value()->GetLargestPossibleRegion() );
          itImageWSX2.GoToBegin();
          itImageWSXX2.GoToBegin();

          for( unsigned int offset=0; offset<j*4; offset++ )
            {
            ++itImageWSX2;
            ++itImageWSXX2;
            }
          for( unsigned int ind=0; ind<4; ind++ )
            {
            m_Data1D( ind )= itImageWSX2.Get();
            ++itImageWSX2;
            }
          itImageWSX2.GoToBegin();
          for( unsigned int offset=0; offset<j; offset++ )
            {
            ++itImageWSX2;
            }
          itImageWSX2.Set( m_Spline1D->DataValue( m_Data1D,
              ( ( x( ( int )m_Dimension-i-1 ) 
              - ( int )x( ( int )m_Dimension-i-1 ) ) ) ) );
          for( unsigned int ind=0; ind<4; ind++ )
            {
            m_Data1D( ind )= itImageWSXX2.Get();
            ++itImageWSXX2;
            }

          itImageWSXX2.GoToBegin();
          for( unsigned int offset=0; offset<j; offset++ )
            {
            ++itImageWSXX2;
            }

          itImageWSXX2.Set( m_Spline1D->DataValue( m_Data1D,
            ( ( x( ( int )m_Dimension-i-1 ) 
            - ( int )x( ( int )m_Dimension-i-1 ) ) ) ) );
          }
        }

      itImageWSX.Set( vD );
      itImageWSXX.Set( vD2 );
      ++itImageWSX;
      ++itImageWSXX;
      }

    itWSX++;
    itWSXX++;
    }

  itWSX  = m_DataWSX->Begin();
  itWSXX = m_DataWSXX->Begin();

  ImageType::IndexType ind;
  for( unsigned int k=0; k<ImageType::ImageDimension; k++ )
    {
    ind[k]=0;
    }

  unsigned int i=0;
  while( itWSX != m_DataWSX->End() )
    {
    d( i ) = itWSX->Value()->GetPixel( ind );
    d2( i ) = itWSXX->Value()->GetPixel( ind );
    itWSX++;
    itWSXX++;
    i++;
    }

  m_Val = m_DataWS->GetPixel( ind );
  return m_Val;
}


bool
SplineND
::Extreme( VectorType & extX, double * extVal )
{
  return m_OptimizerND->Extreme( extX, extVal );
}


bool
SplineND
::Extreme( VectorType & extX, double * extVal, unsigned int n,
  MatrixType & dirs )
{
  return m_OptimizerND->Extreme( extX, extVal, n, dirs );
}


bool
SplineND
::Extreme( VectorType & extX, double * extVal, VectorType & dir )
{
  MatrixType h( m_Dimension, 1 );
  for( unsigned int i=0; i<m_Dimension; i++ )
    {
    h( i, 0 ) = dir( i );
    }
  return m_OptimizerND->Extreme( extX, extVal, 1, h );
}


bool
SplineND
::ExtremeConjGrad( VectorType & extX, double * extVal )
{
  Hessian( extX );

  VectorType eVals( m_Dimension, 0.0 );
  MatrixType eVects( m_Dimension, m_Dimension );
  if( m_OptimizerND->GetSearchForMin() )
    {
    ComputeEigen( m_H, eVects, eVals, /* orderByAbs= */false,
      /* minToMax= */false );
    }
  else
    {
    ComputeEigen( m_H, eVects, eVals, /* orderByAbs= */false,
      /* minToMax= */true );
    }

  return m_OptimizerND->Extreme( extX, extVal, m_Dimension, eVects );
}


void
SplineND
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  os << indent << "Dimension:        " << m_Dimension << std::endl;
  os << indent << "Clip:             " << m_Clip << std::endl;
  os << indent << "XMin:             " << m_XMin << std::endl;
  os << indent << "XMax:             " << m_XMax << std::endl;
  os << indent << "NewData:          " << m_NewData << std::endl;
  os << indent << "Xi:               " << m_Xi << std::endl;
  os << indent << "Val:              " << m_Val << std::endl;
  os << indent << "D:                " << m_D << std::endl;
  os << indent << "H:                " << m_H << std::endl;
  os << indent << "Data:             " << m_Data << std::endl;
  os << indent << "DataWS:           " << m_DataWS << std::endl;
  os << indent << "Data1D:           " << m_Data1D << std::endl;
  os << indent << "DataWSX:          " << m_DataWSX << std::endl;
  os << indent << "DataWSXX:         " << m_DataWSXX << std::endl;
  os << indent << "FuncVal:          " << m_FuncVal << std::endl;
  os << indent << "OptimizerNDVal:   " << m_OptimizerNDVal << std::endl;
  os << indent << "OptimizerNDDeriv: " << m_OptimizerNDDeriv << std::endl;
  os << indent << "OptimizerND:      " << m_OptimizerND << std::endl;
  os << indent << "Spline1D:         " << m_Spline1D << std::endl;
}

} // End namespace tube
