/*=========================================================================

Library:   TubeTK/VTree

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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
#include <cmath>
#include <itkVectorContainer.h>
#include <itkImageRegionIterator.h>

namespace tube
{

class SplineNDValFunc :
public UserFunc< vnl_vector<double>, double >
{
private:

  SplineND * cSpline;
  double cVal;

public:

  SplineNDValFunc( SplineND * newSpline )
    {
    cSpline = newSpline;
    cVal = 0;
    };

  const double & value( const vnl_vector<double> & x )
    {
    cVal = cSpline->value( x );
    return cVal;
    };
};

class SplineNDDerivFunc :
public UserFunc< vnl_vector<double>, vnl_vector<double> >
{
private:

  SplineND * cSpline;
  vnl_vector<double> cDx;

public:

  SplineNDDerivFunc( SplineND * newSpline )
    {
    cSpline = newSpline;
    cDx.set_size( cSpline->nDims() );
    };

  const vnl_vector<double> & value( const vnl_vector<double> & x )
    {
    cDx = cSpline->valueD( x );
    return cDx;
    }
};


SplineND::SplineND()
  {
  m_debug = false;

  cNDims = 0;

  cClip = false;
  cNewData = true;

  cVal = 0;
  cData = NULL;
  cDataWS = NULL;
  cDataWSX = NULL;
  cDataWSXX = NULL;


  cFuncVal = NULL;

  cOptNDVal = new SplineNDValFunc( this );
  cOptNDDeriv = new SplineNDDerivFunc( this );

  cOptND = NULL;
  cSpline1D = NULL;

  this->use( 0, NULL, NULL, NULL );
  }

SplineND::SplineND( unsigned int newNDims,
  UserFunc<IntVectorType, double> * newFuncVal,
  Spline1D * newSpline1D,
  Optimizer1D * newOpt1D )
  {
  m_debug = false;

  cNDims = 0;

  cClip = false;
  cNewData = true;

  cVal = 0;
  cData = NULL;
  cDataWS = NULL;
  cDataWSX = NULL;
  cDataWSXX = NULL;

  cFuncVal = NULL;

  cOptNDVal = new SplineNDValFunc( this );
  cOptNDDeriv = new SplineNDDerivFunc( this );

  cOptND = NULL;
  cSpline1D = NULL;

  this->use( newNDims, newFuncVal, newSpline1D, newOpt1D );
  }

SplineND::~SplineND()
{
  delete cOptNDVal;
  delete cOptNDDeriv;
  if( cOptND != NULL )
    {
    delete cOptND;
    }
}


void SplineND::use( unsigned int newNDims,
  UserFunc< IntVectorType, double > * newFuncVal,
  Spline1D * newSpline1D,
  Optimizer1D * newOpt1D )
  {
  if( m_debug )
    {
    std::cout << "Spline::use()" << std::endl;
    }
  if( cOptND != NULL )
    {
    delete cOptND;
    cOptND = NULL;
    }

  cNDims = newNDims;
  cXMin.set_size( cNDims );
  cXMin.fill( ( int )0 );
  cXMax.set_size( cNDims );
  cXMax.fill( ( int )1 );
  cXi.set_size( cNDims );
  cD.set_size( cNDims );
  cH.set_size( cNDims, cNDims );

  ImageType::SizeType dimSize;
  for( unsigned int i=0; i<dimSize.GetSizeDimension(); i++ )
    {
    dimSize[i] = 1;
    }
  for( unsigned int i=0; i<cNDims; i++ )
    {
    dimSize[i] = 4;
    }

  cData   = ImageType::New();
  ImageType::RegionType region;
  region.SetSize( dimSize );
  cData->SetRegions( region );
  cData->Allocate();

  cDataWS = ImageType::New();
  cDataWS->SetRegions( region );
  cDataWS->Allocate();

  cDataWSX = VectorImageType::New();
  cDataWSXX = VectorImageType::New();
  cDataWSX->Reserve( cNDims );
  cDataWSXX->Reserve( cNDims );

  VectorImageType::Iterator itWSX  = cDataWSX->Begin();
  VectorImageType::Iterator itWSXX = cDataWSXX->Begin();
  while( itWSX != cDataWSX->End() )
    {
    itWSX->Value() = ImageType::New();
    ImageType::RegionType region;

    region.SetSize( dimSize );
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

  cData1D.set_size( 4 );

  cVal = 0;
  cFuncVal = newFuncVal;
  cSpline1D = newSpline1D;

  if( newOpt1D != NULL )
    {
    cOptND = new OptimizerND( cNDims, cOptNDVal, cOptNDDeriv, newOpt1D );
    }

  cNewData = true;
  }

//
//
//
bool SplineND::clipEdge( void )
  {
  return cClip;
  }

void SplineND::clipEdge( bool newClip )
  {
  cClip = newClip;

  }

//
//
//
const SplineND::IntVectorType & SplineND::xMin( void )
  {
  return cXMin;
  }

void SplineND::xMin( const IntVectorType & newXMin )
  {
  VectorType t( cNDims );
  for( unsigned int i=0; i<cNDims; i++ )
    {
    t( i ) = newXMin( i );
    }
  cOptND->xMin( t );
  cXMin = newXMin;
  }


const SplineND::IntVectorType & SplineND::xMax( void )
  {
  return cXMax;
  }

void SplineND::xMax( const IntVectorType & newXMax )
  {
  VectorType t( cNDims );
  for( unsigned int i=0; i<cNDims; i++ )
    {
    t( i ) = newXMax( i );
    }
  cOptND->xMax( t );
  cXMax = newXMax;
  }


//
//
//
bool SplineND::newData( void )
  {
  return cNewData;
  }

void SplineND::newData( bool newNewData )
  {
  cNewData = newNewData;
  }

//
//
//
void SplineND::cGetData( const VectorType & x )
  {
  if( m_debug )
    {
    std::cout << "SplineND: cGetData: " << x << std::endl;
    }

  bool eql = true;
  if( !cNewData )
    {
    for( unsigned int i=0; i<cNDims; i++ )
      {
      if( cXi( i ) != ( int )x( i ) )
        {
        eql = false;
        break;
        }
      }
    }
  if( !eql || cNewData )
    {
    if( cNewData )
      {
      cNewData = false;
      for( unsigned int i=0; i<cNDims; i++ )
        {
        cXi( i ) = ( int )x( i );
        }

      itk::ImageRegionIterator<ImageType> it( cData,
        cData->GetLargestPossibleRegion() );
      it.GoToBegin();

      IntVectorType p( cNDims );
      IntVectorType xiOffset( cNDims, -1 );
      bool done = false;
      while( !done )
        {
        p = cXi + xiOffset;
        for( unsigned int i=0; i<cNDims; i++ )
          {
          if( p( i ) < cXMin( i ) )
            {
            if( cClip )
              {
              p( i ) = cXMin( i );
              }
            else
              {
              p( i ) = cXMin( i ) + ( cXMin( i ) - p( i ) );
              if( p( i ) > cXMax( i ) )
                {
                p( i ) = cXMax( i );
                }
              }
            }
          else if( p( i ) > cXMax( i ) )
            {
            if( cClip )
              {
              p( i ) = cXMax( i );
              }
            else
              {
              p( i ) = cXMax( i ) - ( p( i ) - cXMax( i ) );
              if( p( i ) < cXMin( i ) )
                {
                p( i ) = cXMin( i );
                }
              }
            }
          }
        it.Set( cFuncVal->value( p ) );
        //if( m_debug )
          //{
          //std::cout << "cData: " << p( 0 );
          //for( unsigned int i=1; i<cNDims; i++ )
            //{
            //std::cout << ", " << p( i );
            //}
          //std::cout << " ( " << xiOffset( 0 );
          //for( unsigned int i=1; i<cNDims; i++ )
            //{
            //std::cout << ", " << xiOffset( i );
            //}
          //std::cout << " ) = " << it.Get() << std::endl;
          //}
        ++it;
        unsigned int dim = 0;
        while( !done && dim<cNDims && ( ++xiOffset( dim ) )>2 )
          {
          xiOffset( dim++ ) = -1;
          if( dim >= cNDims )
            {
            done = true;
            }
          }
        }
      }
    else
      {
      IntVectorType p( cNDims );
      IntVectorType pOld( cNDims );
      IntVectorType xiOffset( cNDims );
      IntVectorType shift( cNDims, 100.0 );

      for( unsigned int i=0; i<cNDims; i++ )
        {
        shift[i] = ( int )x( i ) - cXi( i );
        }

      for( unsigned int i=0; i<cNDims; i++ )
        {
        cXi( i ) = ( int )x( i );
        }

      itk::ImageRegionIterator<ImageType> it( cDataWS,
        cDataWS->GetLargestPossibleRegion() );
      it.GoToBegin();

      xiOffset = -1;
      bool done = false;
      while( !done )
        {
        p = cXi + xiOffset;
        pOld = xiOffset + shift + 1;
        //if( m_debug )
          //{
          //std::cout << "newDataPoint=" << xiOffset+1
            //<< " : oldDataPoint=" << pOld << std::endl;
          //}
        for( unsigned int i=0; i<cNDims; i++ )
          {
          if( p( i ) < cXMin( i ) )
            {
            if( cClip )
              {
              p( i ) = cXMin( i );
              }
            else
              {
              p( i ) = cXMin( i ) + ( cXMin( i ) - p( i ) );
              if( p( i ) > cXMax( i ) )
                {
                p( i ) = cXMax( i );
                }
              }
            }
          else if( p( i ) > cXMax( i ) )
            {
            if( cClip )
              {
              p( i ) = cXMax( i );
              }
            else
              {
              p( i ) = cXMax( i ) - ( p( i ) - cXMax( i ) );
              if( p( i ) < cXMin( i ) )
                {
                p( i ) = cXMin( i );
                }
              }
            }
          }
        bool reuse = true;
        for( unsigned int j=0; j<cNDims; j++ )
          {
          if( pOld[j] < 0 || pOld[j] > 3 )
            {
            reuse = false;
            break;
            }
          }
        if( reuse )
          {
          if( m_debug )
            {
            std::cout << "reusing" << std::endl;
            }
          ImageType::IndexType indx;
          indx.Fill( 0 );
          for( unsigned int i=0; i<cNDims; i++ )
            {
            indx[i] = pOld[i];
            }
          it.Set( cData->GetPixel( indx ) );
          }
        else
          {
          it.Set( cFuncVal->value( p ) );
          //if( m_debug )
            //{
            //std::cout << "adding cData: " << p( 0 );
            //for( unsigned int i=1; i<cNDims; i++ )
              //{
              //std::cout << ", " << p( i );
              //}
            //std::cout << " ( " << xiOffset( 0 );
            //for( unsigned int i=1; i<cNDims; i++ )
              //{
              //std::cout << ", " << xiOffset( i );
              //}
            //std::cout << " ) = " << it.Get() << std::endl;
            //}
          }
        ++it;
        unsigned int dim = 0;
        while( !done && dim<cNDims && ( ++xiOffset( dim ) )>2 )
          {
          xiOffset( dim++ ) = -1;
          if( dim >= cNDims )
            {
            done = true;
            }
          }
        }

      if( m_debug )
        {
        std::cout << "cData = " << std::endl;
        }
      itk::ImageRegionIterator<ImageType> itDest( cData,
        cData->GetLargestPossibleRegion() );
      it.GoToBegin();
      itDest.GoToBegin();
      while( !it.IsAtEnd() )
        {
        itDest.Set( it.Get() );
        //if( m_debug )
          //{
          //std::cout << " " << it.GetIndex()
            //<< " = " << it.Get() << std::endl;
          //}
        ++it;
        ++itDest;
        }
      if( m_debug )
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
  if( m_debug )
    {
    std::cout << "SplineND::value()" << std::endl;
    }

  this->cGetData( x );

  itk::ImageRegionIterator<ImageType> itData( cData,
    cData->GetLargestPossibleRegion() );

  itk::ImageRegionIterator<ImageType> itDataWS( cDataWS,
    cDataWS->GetLargestPossibleRegion() );

  itData.GoToBegin();
  itDataWS.GoToBegin();
  while( !itData.IsAtEnd() )
    {
    itDataWS.Set( itData.Get() );
    ++itDataWS;
    ++itData;
    }

  itk::ImageRegionIterator<ImageType> itDataWSColumn( cDataWS,
    cDataWS->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> itDataWSDest( cDataWS,
    cDataWS->GetLargestPossibleRegion() );

  for( int i=( int )cNDims-1; i>=0; i-- )
    {
    unsigned int k = ( unsigned int )pow( ( float )4, ( int )i );
    itDataWSColumn.GoToBegin();
    itDataWSDest.GoToBegin();
    for( unsigned int j=0; j<k; j++ )
      {
      for( unsigned int ind=0; ind<4; ind++ )
        {
        cData1D( ind )= itDataWSColumn.Get();
        ++itDataWSColumn;
        }

      itDataWSDest.Set( cSpline1D->dataValue( cData1D,
          ( ( x( ( int )cNDims-i-1 )-( int )x( ( int )cNDims-i-1 ) ) ) ) );
      ++itDataWSDest;
      }
    }

  itDataWS.GoToBegin();
  if( m_debug )
    {
    std::cout << "SplineND : value : value at " << x << " = "
      << itDataWS.Get() << std::endl;
    }
  cVal = itDataWS.Get();
  return cVal;
  }

//
//
//
double SplineND::valueD( const VectorType & x, IntVectorType & dx )
  {
  this->cGetData( x );

  itk::ImageRegionIterator<ImageType> itData( cData,
    cData->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> itDataWS( cDataWS,
    cDataWS->GetLargestPossibleRegion() );

  itData.GoToBegin();
  itDataWS.GoToBegin();
  while( !itData.IsAtEnd() )
    {
    itDataWS.Set( itData.Get() );
    ++itDataWS;
    ++itData;
    }

  itk::ImageRegionIterator<ImageType> itDataWSColumn( cDataWS,
    cDataWS->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> itDataWSDest( cDataWS,
    cDataWS->GetLargestPossibleRegion() );

  for( int i=( int )cNDims-1; i>=0; i-- )
    {
    itDataWSColumn.GoToBegin();
    itDataWSDest.GoToBegin();
    unsigned int k = ( unsigned int )pow( ( float )4, ( int )i );
    switch( dx( ( int )cNDims-i-1 ) )
      {
      default:
      case 0:
        for( unsigned int j=0; j<k; j++ )
          {
          for( unsigned int ind=0; ind<4; ind++ )
            {
            cData1D( ind )= itDataWSColumn.Get();
            ++itDataWSColumn;
            }

          itDataWSDest.Set( cSpline1D->dataValue( cData1D,
              ( ( x( ( int )cNDims-i-1 )-( int )x( ( int )cNDims-i-1 ) ) ) ) );
          ++itDataWSDest;
          }
        break;
      case 1:
        for( unsigned int j=0; j<k; j++ )
          {
          for( unsigned int ind=0; ind<4; ind++ )
            {
            cData1D( ind )= itDataWSColumn.Get();
            ++itDataWSColumn;
            }

          itDataWSDest.Set( cSpline1D->dataValueD( cData1D,
              ( ( x( ( int )cNDims-i-1 )-( int )x( ( int )cNDims-i-1 ) ) ) ) );
          ++itDataWSDest;
          }
        break;
      case 2:
        for( unsigned int j=0; j<k; j++ )
          {
          for( unsigned int ind=0; ind<4; ind++ )
            {
            cData1D( ind )= itDataWSColumn.Get();
            ++itDataWSColumn;
            }

          itDataWSDest.Set( cSpline1D->dataValueD2( cData1D,
              ( ( x( ( int )cNDims-i-1 )-( int )x( ( int )cNDims-i-1 ) ) ) ) );
          ++itDataWSDest;
          }
        break;
      }
    }

  itDataWS.GoToBegin();
  cVal = itDataWS.Get();
  return cVal;
  }


SplineND::VectorType & SplineND::valueD( const VectorType & x )
  {
  IntVectorType dx( cNDims, 0 );

  for( unsigned int i=0; i<cNDims; i++ )
    {
    dx( i ) = 1;
    cD( i ) = valueD( x, dx );
    dx( i ) = 0;
    }

  return cD;
  }

//
//
//
SplineND::MatrixType & SplineND::hessian( const VectorType & x )
  {

  IntVectorType dx( cNDims, 0 );

  VectorType d( cNDims );
  VectorType d2( cNDims );
  valueVDD2( x, d, d2 );
  for( unsigned int i=0; i<cNDims; i++ )
    {
    cH.put( i, i, d2( i ) );
    cD( i ) = d( i );
    }

  for( unsigned int i=0; i<cNDims; i++ )
    {
    for( unsigned int j=i+1; j<cNDims; j++ )
      {
      dx( i ) = 1;
      dx( j ) = 1;
      cH.put( i, j, valueD( x, dx ) );
      cH.put( j, i, cH.get( i, j ) );
      dx( i ) = 0;
      dx( j ) = 0;
      }
    }
  return cH;
  }

//
//
//
double SplineND::valueJet( const VectorType & x,
  VectorType & d, MatrixType & h )
  {
  IntVectorType dx( cNDims, 0 );

  VectorType lD( cNDims );
  VectorType lD2( cNDims );
  double v = valueVDD2( x, lD, lD2 );
  MatrixType tempMatrix( cNDims, cNDims );

  if( m_debug )
    {
    std::cout << "SplineND::valueJet() v= " << v << std::endl;
    }
  for( unsigned int i=0; i< cNDims; i++ )
    {
    cH.put( i, i, lD2( i ) );
    cD( i ) = lD( i );
    }



  for( unsigned int i=0; i<cNDims; i++ )
    {
    for( unsigned int j=i+1; j<cNDims; j++ )
      {
      dx( i ) = 1;
      dx( j ) = 1;
      cH.put( i, j, valueD( x, dx ) );
      cH.put( j, i, cH.get( i, j ) );
      dx( i ) = 0;
      dx( j ) = 0;
      }
    }

  if( m_debug )
    {
    std::cout << "SplineND : valueJet : value at " << x
      << " = " << v << std::endl;
    }

  for( unsigned int i=0; i<cNDims; i++ )
    {
    d( i ) = cD( i );
    }

  h = cH;

  return v;
  }

//
//
double SplineND::valueVDD2( const VectorType & x,
  VectorType & d, VectorType & d2 )
  {
  this->cGetData( x );

  VectorImageType::Iterator itWSX  = cDataWSX->Begin();
  VectorImageType::Iterator itWSXX = cDataWSXX->Begin();

  itk::ImageRegionIterator<ImageType> itData( cData,
    cData->GetLargestPossibleRegion() );
  itk::ImageRegionIterator<ImageType> itDataWS( cDataWS,
    cDataWS->GetLargestPossibleRegion() );

  /** Could optimize that */
  itDataWS.GoToBegin();
  itData.GoToBegin();
  while( !itData.IsAtEnd() )
    {
    itDataWS.Set( itData.Get() );
    ++itDataWS;
    ++itData;
    }

  while( itWSX != cDataWSX->End() )
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

  itWSX  = cDataWSX->Begin();
  itWSXX = cDataWSXX->Begin();
  double vD, vD2;

  for( int i=( int )cNDims-1; i>=0; i-- )
    {
    unsigned int k = ( unsigned int )pow( ( float )4, ( int )i );

    itk::ImageRegionIterator<ImageType> itImageWSX( itWSX->Value(),
      itWSX->Value()->GetLargestPossibleRegion() );
    itk::ImageRegionIterator<ImageType> itImageWSXX( itWSXX->Value(),
      itWSXX->Value()->GetLargestPossibleRegion() );
    itImageWSX.GoToBegin();
    itImageWSXX.GoToBegin();

    itk::ImageRegionIterator<ImageType> itDataWSX( cDataWS,
      cDataWS->GetLargestPossibleRegion() );
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
        cData1D( ind )= itDataWSX.Get();
        ++itDataWSX;
        }

      itDataWSX.GoToBegin();
      for( unsigned int offset=0; offset<j; offset++ )
        {
        ++itDataWSX;
        }

      itDataWSX.Set( cSpline1D->dataValueJet( cData1D,
        ( ( x( ( int )cNDims-i-1 )-( int )x( ( int )cNDims-i-1 ) ) ), &vD, &vD2 ) );

      for( unsigned int l=0; l<cNDims; l++ )
        {

        if( ( int )cNDims-i != ( int )l )
          {
          VectorImageType::Iterator itWSX2  = cDataWSX->Begin();
          VectorImageType::Iterator itWSXX2 = cDataWSXX->Begin();

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
            cData1D( ind )= itImageWSX2.Get();
            ++itImageWSX2;
            }
          itImageWSX2.GoToBegin();
          for( unsigned int offset=0; offset<j; offset++ )
            {
            ++itImageWSX2;
            }
          itImageWSX2.Set( cSpline1D->dataValue( cData1D,
              ( ( x( ( int )cNDims-i-1 )-( int )x( ( int )cNDims-i-1 ) ) ) ) );
          for( unsigned int ind=0; ind<4; ind++ )
            {
            cData1D( ind )= itImageWSXX2.Get();
            ++itImageWSXX2;
            }

          itImageWSXX2.GoToBegin();
          for( unsigned int offset=0; offset<j; offset++ )
            {
            ++itImageWSXX2;
            }

          itImageWSXX2.Set( cSpline1D->dataValue( cData1D,
            ( ( x( ( int )cNDims-i-1 )-( int )x( ( int )cNDims-i-1 ) ) ) ) );
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

  itWSX  = cDataWSX->Begin();
  itWSXX = cDataWSXX->Begin();

  ImageType::IndexType ind;
  for( unsigned int k=0; k<ImageType::ImageDimension; k++ )
    {
    ind[k]=0;
    }

  unsigned int i=0;
  while( itWSX != cDataWSX->End() )
    {
    d( i ) = itWSX->Value()->GetPixel( ind );
    d2( i ) = itWSXX->Value()->GetPixel( ind );
    itWSX++;
    itWSXX++;
    i++;
    }

  cVal = cDataWS->GetPixel( ind );
  return cVal;
  }

//
//
//
OptimizerND * SplineND::optimizerND( void )
  {
  return cOptND;
  }

//
//
//
bool SplineND::extreme( VectorType & extX, double * extVal )
  {
  return cOptND->extreme( extX, extVal );
  }

bool SplineND::extreme( VectorType & extX, double *extVal,
  unsigned int n, MatrixType &dirs )
  {
  return cOptND->extreme( extX, extVal, n, dirs );
  }

bool SplineND::extreme( VectorType & extX, double *extVal, VectorType &dir )
  {
  MatrixType h( cNDims, 1 );
  for( unsigned int i=0; i<cNDims; i++ )
    {
    h( i,0 ) = dir( i );
    }
  return cOptND->extreme( extX, extVal, 1, h );
  }

bool SplineND::extremeConjGrad( VectorType & extX, double * extVal )
  {
  hessian( extX );

  //  VectorType eVals( cNDims );
  VectorType eVals( cNDims, NULL );
  MatrixType eVects( cNDims, cNDims );
  Eigen( cH, eVects, eVals, false );

  return cOptND->extreme( extX, extVal, cNDims, eVects );
  }


}; // namespace
