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
#include "tubeSplineND.h"
#include "tubeUserFunc.h"
#include <cmath>
#include <itkVectorContainer.h>
#include <itkImageRegionIterator.h>

namespace tube
{

class SplineNDValFunc : public UserFunc< vnl_vector<double>, double >
{
public:
  SplineNDValFunc( SplineND * newSpline )
    {
    m_Spline = newSpline;
    m_Val = 0;
    }

  const double & value( const vnl_vector<double> & x )
    {
    m_Val = m_Spline->value( x );
    return m_Val;
    }

private:
  SplineND * m_Spline;
  double     m_Val;

}; // End class SplineNDValFunc

class SplineNDDerivFunc : public UserFunc< vnl_vector<double>, vnl_vector<double> >
{
public:
  SplineNDDerivFunc( SplineND * newSpline )
    {
    m_Spline = newSpline;
    m_Dx.set_size( m_Spline->nDims() );
    }

  const vnl_vector<double> & value( const vnl_vector<double> & x )
    {
    m_Dx = m_Spline->valueD( x );
    return m_Dx;
    }

private:
  SplineND * m_Spline;
  vnl_vector<double> m_Dx;

}; // End class SplineNDDerivFunc


SplineND::SplineND( void )
{
  m_Debug = false;

  m_NDims = 0;

  m_Clip = false;
  m_NewData = true;

  m_Val = 0;
  m_Data = NULL;
  m_DataWS = NULL;
  m_DataWSX = NULL;
  m_DataWSXX = NULL;


  m_FuncVal = NULL;

  m_OptNDVal = new SplineNDValFunc( this );
  m_OptNDDeriv = new SplineNDDerivFunc( this );

  m_OptND = NULL;
  m_Spline1D = NULL;

  this->use( 0, NULL, NULL, NULL );
}

SplineND::SplineND( unsigned int newNDims,
  UserFunc<IntVectorType, double> * newFunm_Val,
  Spline1D * newSpline1D,
  Optimizer1D * newOpt1D )
{
  m_Debug = false;

  m_NDims = 0;

  m_Clip = false;
  m_NewData = true;

  m_Val = 0;
  m_Data = NULL;
  m_DataWS = NULL;
  m_DataWSX = NULL;
  m_DataWSXX = NULL;

  m_FuncVal = NULL;

  m_OptNDVal = new SplineNDValFunc( this );
  m_OptNDDeriv = new SplineNDDerivFunc( this );

  m_OptND = NULL;
  m_Spline1D = NULL;

  this->use( newNDims, newFunm_Val, newSpline1D, newOpt1D );
}

SplineND::~SplineND( void )
{
  delete m_OptNDVal;
  delete m_OptNDDeriv;
  if( m_OptND != NULL )
    {
    delete m_OptND;
    }
}


void SplineND::use( unsigned int newNDims,
  UserFunc< IntVectorType, double > * newFunm_Val,
  Spline1D * newSpline1D,
  Optimizer1D * newOpt1D )
{
  if( m_Debug )
    {
    std::cout << "Spline::use()" << std::endl;
    }
  if( m_OptND != NULL )
    {
    delete m_OptND;
    m_OptND = NULL;
    }

  m_NDims = newNDims;
  m_XMin.set_size( m_NDims );
  m_XMin.fill( ( int )0 );
  m_XMax.set_size( m_NDims );
  m_XMax.fill( ( int )1 );
  m_Xi.set_size( m_NDims );
  m_D.set_size( m_NDims );
  m_H.set_size( m_NDims, m_NDims );

  ImageType::SizeType dimSize;
  for( unsigned int i=0; i<dimSize.GetSizeDimension(); i++ )
    {
    dimSize[i] = 1;
    }
  for( unsigned int i=0; i<m_NDims; i++ )
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
  m_DataWSX->Reserve( m_NDims );
  m_DataWSXX->Reserve( m_NDims );

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
  m_FuncVal = newFunm_Val;
  m_Spline1D = newSpline1D;

  if( newOpt1D != NULL )
    {
    m_OptND = new OptimizerND( m_NDims, m_OptNDVal, m_OptNDDeriv, newOpt1D );
    }

  m_NewData = true;
}

//
//
//
bool SplineND::clipEdge( void )
{
  return m_Clip;
}

void SplineND::clipEdge( bool newClip )
{
  m_Clip = newClip;
}

//
//
//
const SplineND::IntVectorType & SplineND::xMin( void )
{
  return m_XMin;
}

void SplineND::xMin( const IntVectorType & newXMin )
{
  VectorType t( m_NDims );
  for( unsigned int i=0; i<m_NDims; i++ )
    {
    t( i ) = newXMin( i );
    }
  m_OptND->xMin( t );
  m_XMin = newXMin;
}


const SplineND::IntVectorType & SplineND::xMax( void )
{
  return m_XMax;
}

void SplineND::xMax( const IntVectorType & newXMax )
{
  VectorType t( m_NDims );
  for( unsigned int i=0; i<m_NDims; i++ )
    {
    t( i ) = newXMax( i );
    }
  m_OptND->xMax( t );
  m_XMax = newXMax;
}


//
//
//
bool SplineND::newData( void )
{
  return m_NewData;
}

void SplineND::newData( bool newNewData )
{
  m_NewData = newNewData;
}

//
//
//
void SplineND::m_GetData( const VectorType & x )
{
  if( m_Debug )
    {
    std::cout << "SplineND: m_GetData: " << x << std::endl;
    }

  bool eql = true;
  if( !m_NewData )
    {
    for( unsigned int i=0; i<m_NDims; i++ )
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
      for( unsigned int i=0; i<m_NDims; i++ )
        {
        m_Xi( i ) = ( int )x( i );
        }

      itk::ImageRegionIterator<ImageType> it( m_Data,
        m_Data->GetLargestPossibleRegion() );
      it.GoToBegin();

      IntVectorType p( m_NDims );
      IntVectorType xiOffset( m_NDims, -1 );
      bool done = false;
      while( !done )
        {
        p = m_Xi + xiOffset;
        for( unsigned int i=0; i<m_NDims; i++ )
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
        it.Set( m_FuncVal->value( p ) );
        if( m_Debug )
          {
          std::cout << "m_Data: " << p( 0 );
          for( unsigned int i=1; i<m_NDims; i++ )
            {
            std::cout << ", " << p( i );
            }
          std::cout << " ( " << xiOffset( 0 );
          for( unsigned int i=1; i<m_NDims; i++ )
            {
            std::cout << ", " << xiOffset( i );
            }
          std::cout << " ) = " << it.Get() << std::endl;
          }
        ++it;
        unsigned int dim = 0;
        while( !done && dim<m_NDims && ( ++xiOffset( dim ) )>2 )
          {
          xiOffset( dim++ ) = -1;
          if( dim >= m_NDims )
            {
            done = true;
            }
          }
        }
      }
    else
      {
      IntVectorType p( m_NDims );
      IntVectorType pOld( m_NDims );
      IntVectorType xiOffset( m_NDims );
      IntVectorType shift( m_NDims, 100.0 );

      for( unsigned int i=0; i<m_NDims; i++ )
        {
        shift[i] = ( int )x( i ) - m_Xi( i );
        }

      for( unsigned int i=0; i<m_NDims; i++ )
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
        if( m_Debug )
          {
          std::cout << "newDataPoint=" << xiOffset+1
            << " : oldDataPoint=" << pOld << std::endl;
          }
        for( unsigned int i=0; i<m_NDims; i++ )
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
        for( unsigned int j=0; j<m_NDims; j++ )
          {
          if( pOld[j] < 0 || pOld[j] > 3 )
            {
            reuse = false;
            break;
            }
          }
        if( reuse )
          {
          if( m_Debug )
            {
            std::cout << "reusing" << std::endl;
            }
          ImageType::IndexType indx;
          indx.Fill( 0 );
          for( unsigned int i=0; i<m_NDims; i++ )
            {
            indx[i] = pOld[i];
            }
          it.Set( m_Data->GetPixel( indx ) );
          }
        else
          {
          it.Set( m_FuncVal->value( p ) );
          if( m_Debug )
            {
            std::cout << "adding m_Data: " << p( 0 );
            for( unsigned int i=1; i<m_NDims; i++ )
              {
              std::cout << ", " << p( i );
              }
            std::cout << " ( " << xiOffset( 0 );
            for( unsigned int i=1; i<m_NDims; i++ )
              {
              std::cout << ", " << xiOffset( i );
              }
            std::cout << " ) = " << it.Get() << std::endl;
            }
          }
        ++it;
        unsigned int dim = 0;
        while( !done && dim<m_NDims && ( ++xiOffset( dim ) )>2 )
          {
          xiOffset( dim++ ) = -1;
          if( dim >= m_NDims )
            {
            done = true;
            }
          }
        }

      if( m_Debug )
        {
        std::cout << "m_Data = " << std::endl;
        }
      itk::ImageRegionIterator<ImageType> itDest( m_Data,
        m_Data->GetLargestPossibleRegion() );
      it.GoToBegin();
      itDest.GoToBegin();
      while( !it.IsAtEnd() )
        {
        itDest.Set( it.Get() );
        if( m_Debug )
          {
          std::cout << " " << it.GetIndex()
            << " = " << it.Get() << std::endl;
          }
        ++it;
        ++itDest;
        }
      if( m_Debug )
        {
        std::cout << "...done" << std::endl;
        }
      }
    }
}


//
//
//
const double & SplineND::value( const VectorType & x )
{
  if( m_Debug )
    {
    std::cout << "SplineND::value()" << std::endl;
    }

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

  for( int i=( int )m_NDims-1; i>=0; i-- )
    {
    unsigned int k = ( unsigned int )pow( ( float )4, ( int )i );
    itDataWSColumn.GoToBegin();
    itDataWSDest.GoToBegin();
    for( unsigned int j=0; j<k; j++ )
      {
      for( unsigned int ind=0; ind<4; ind++ )
        {
        m_Data1D( ind )= itDataWSColumn.Get();
        ++itDataWSColumn;
        }

      itDataWSDest.Set( m_Spline1D->dataValue( m_Data1D,
          ( ( x( ( int )m_NDims-i-1 )-( int )x( ( int )m_NDims-i-1 ) ) ) ) );
      ++itDataWSDest;
      }
    }

  itDataWS.GoToBegin();
  if( m_Debug )
    {
    std::cout << "SplineND : value : value at " << x << " = "
      << itDataWS.Get() << std::endl;
    }
  m_Val = itDataWS.Get();
  return m_Val;
}

//
//
//
double SplineND::valueD( const VectorType & x, IntVectorType & dx )
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

  for( int i=( int )m_NDims-1; i>=0; i-- )
    {
    itDataWSColumn.GoToBegin();
    itDataWSDest.GoToBegin();
    unsigned int k = ( unsigned int )pow( ( float )4, ( int )i );
    switch( dx( ( int )m_NDims-i-1 ) )
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

          itDataWSDest.Set( m_Spline1D->dataValue( m_Data1D,
              ( ( x( ( int )m_NDims-i-1 )-( int )x( ( int )m_NDims-i-1 ) ) ) ) );
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

          itDataWSDest.Set( m_Spline1D->dataValueD( m_Data1D,
              ( ( x( ( int )m_NDims-i-1 )-( int )x( ( int )m_NDims-i-1 ) ) ) ) );
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

          itDataWSDest.Set( m_Spline1D->dataValueD2( m_Data1D,
              ( ( x( ( int )m_NDims-i-1 )-( int )x( ( int )m_NDims-i-1 ) ) ) ) );
          ++itDataWSDest;
          }
        break;
      }
    }

  itDataWS.GoToBegin();
  m_Val = itDataWS.Get();
  return m_Val;
}


SplineND::VectorType & SplineND::valueD( const VectorType & x )
{
  IntVectorType dx( m_NDims, 0 );

  for( unsigned int i=0; i<m_NDims; i++ )
    {
    dx( i ) = 1;
    m_D( i ) = valueD( x, dx );
    dx( i ) = 0;
    }

  return m_D;
}

//
//
//
SplineND::MatrixType & SplineND::hessian( const VectorType & x )
{

  IntVectorType dx( m_NDims, 0 );

  VectorType d( m_NDims );
  VectorType d2( m_NDims );
  valueVDD2( x, d, d2 );
  for( unsigned int i=0; i<m_NDims; i++ )
    {
    m_H.put( i, i, d2( i ) );
    m_D( i ) = d( i );
    }

  for( unsigned int i=0; i<m_NDims; i++ )
    {
    for( unsigned int j=i+1; j<m_NDims; j++ )
      {
      dx( i ) = 1;
      dx( j ) = 1;
      m_H.put( i, j, valueD( x, dx ) );
      m_H.put( j, i, m_H.get( i, j ) );
      dx( i ) = 0;
      dx( j ) = 0;
      }
    }
  return m_H;
}

//
//
//
double SplineND::valueJet( const VectorType & x,
  VectorType & d, MatrixType & h )
{
  IntVectorType dx( m_NDims, 0 );

  VectorType lD( m_NDims );
  VectorType lD2( m_NDims );
  double v = valueVDD2( x, lD, lD2 );
  MatrixType tempMatrix( m_NDims, m_NDims );

  if( m_Debug )
    {
    std::cout << "SplineND::valueJet() v= " << v << std::endl;
    }
  for( unsigned int i=0; i< m_NDims; i++ )
    {
    m_H.put( i, i, lD2( i ) );
    m_D( i ) = lD( i );
    }

  for( unsigned int i=0; i<m_NDims; i++ )
    {
    for( unsigned int j=i+1; j<m_NDims; j++ )
      {
      dx( i ) = 1;
      dx( j ) = 1;
      m_H.put( i, j, valueD( x, dx ) );
      m_H.put( j, i, m_H.get( i, j ) );
      dx( i ) = 0;
      dx( j ) = 0;
      }
    }

  if( m_Debug )
    {
    std::cout << "SplineND : valueJet : value at " << x
      << " = " << v << std::endl;
    }

  for( unsigned int i=0; i<m_NDims; i++ )
    {
    d( i ) = m_D( i );
    }

  h = m_H;

  return v;
}

//
//
double SplineND::valueVDD2( const VectorType & x,
  VectorType & d, VectorType & d2 )
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

  for( int i=( int )m_NDims-1; i>=0; i-- )
    {
    unsigned int k = ( unsigned int )pow( ( float )4, ( int )i );

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

      itDataWSX.Set( m_Spline1D->dataValueJet( m_Data1D,
        ( ( x( ( int )m_NDims-i-1 )-( int )x( ( int )m_NDims-i-1 ) ) ), &vD, &vD2 ) );

      for( unsigned int l=0; l<m_NDims; l++ )
        {

        if( ( int )m_NDims-i != ( int )l )
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
          itImageWSX2.Set( m_Spline1D->dataValue( m_Data1D,
              ( ( x( ( int )m_NDims-i-1 )-( int )x( ( int )m_NDims-i-1 ) ) ) ) );
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

          itImageWSXX2.Set( m_Spline1D->dataValue( m_Data1D,
            ( ( x( ( int )m_NDims-i-1 )-( int )x( ( int )m_NDims-i-1 ) ) ) ) );
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

//
//
//
OptimizerND * SplineND::optimizerND( void )
{
  return m_OptND;
}

//
//
//
bool SplineND::extreme( VectorType & extX, double * extVal )
{
  return m_OptND->extreme( extX, extVal );
}

bool SplineND::extreme( VectorType & extX, double *extVal,
  unsigned int n, MatrixType &dirs )
{
  return m_OptND->extreme( extX, extVal, n, dirs );
}

bool SplineND::extreme( VectorType & extX, double *extVal, VectorType &dir )
{
  MatrixType h( m_NDims, 1 );
  for( unsigned int i=0; i<m_NDims; i++ )
    {
    h( i,0 ) = dir( i );
    }
  return m_OptND->extreme( extX, extVal, 1, h );
}

bool SplineND::extremeConjGrad( VectorType & extX, double * extVal )
{
  hessian( extX );

  VectorType eVals( m_NDims, 0.0 );
  MatrixType eVects( m_NDims, m_NDims );
  ::tube::ComputeEigen( m_H, eVects, eVals, false );

  return m_OptND->extreme( extX, extVal, m_NDims, eVects );
}

} // End namespace tube
