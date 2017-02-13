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

#include "tubeBrentOptimizer1D.h"
#include "tubeSplineApproximation1D.h"

#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

bool MY_DEBUG = false;

class MySA1DFunc : public tube::UserFunction< int, double >
{
private:
  double cVal;

public:
  MySA1DFunc( void )
    {
    cVal = 0;
    }
  const double & Value( const int & x )
    {
    cVal = std::sin( ( double )x );
    if( MY_DEBUG )
      {
      std::cout << "   x = " << x << " : v = " << cVal << std::endl;
      }
    return cVal;
    }

}; // End class MySA1DFunc

class MySA1DFuncV : public tube::UserFunction< double, double >
{
private:
  double cVal;

public:
  MySA1DFuncV( void )
    {
    cVal = 0;
    }
  const double & Value( const double & x )
    {
    cVal = std::sin( ( double )x );
    if( MY_DEBUG )
      {
      std::cout << "   x = " << x << " : v = " << cVal << std::endl;
      }
    return cVal;
    }

}; // End class MySA1DFuncV

class MySA1DFuncD : public tube::UserFunction< double, double >
{
private:
  double cDeriv;

public:
  MySA1DFuncD( void )
    {
    cDeriv = 0;
    }
  const double & Value( const double & x )
    {
    cDeriv = std::cos( ( double )x );
    if( MY_DEBUG )
      {
      std::cout << "   x = " << x << " : dx = " << cDeriv << std::endl;
      }
    return cDeriv;
    }

}; // End class MySA1DFuncD

int tubeSplineApproximation1DTest( int argc, char * argv[] )
{
  if( argc != 2 )
    {
    std::cout << "usage: run <outImFile>" << std::endl;
    return EXIT_FAILURE;
    }

  double epsilon = 0.000001;

  MySA1DFunc * myFunc = new MySA1DFunc();
  MySA1DFuncV * myFuncV = new MySA1DFuncV();
  MySA1DFuncD * myFuncD = new MySA1DFuncD();

  tube::BrentOptimizer1D * opt = new tube::BrentOptimizer1D();
  opt->SetEpsilon( epsilon );
  opt->SetSearchForMin( true );
  opt->SetXStep( 0.01 );
  opt->SetTolerance( 0.0000001 );

  tube::SplineApproximation1D spline( myFunc, opt );

  int returnStatus = EXIT_SUCCESS;

  spline.SetClip( true );

  spline.SetXMin( -3 );
  if( spline.GetXMin() != -3 )
    {
    std::cout << "xMin should be -3 and not " << spline.GetXMin() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  spline.SetXMax( 6 );
  if( spline.GetXMax() != 6 )
    {
    std::cout << "xMax should be 6 and not " << spline.GetXMax() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  typedef itk::Image< float, 3 >  ImageType;

  ImageType::Pointer im = ImageType::New();
  ImageType::RegionType imRegion;
  ImageType::SizeType imSize;
  imSize[0] = 40;
  imSize[1] = 40;
  imSize[2] = 20;
  imRegion.SetSize( imSize );
  ImageType::IndexType index0;
  index0[0] = -10;
  index0[1] = -10;
  index0[2] = -5;
  imRegion.SetIndex( index0 );
  im->SetRegions( imRegion );
  ImageType::SpacingType imSpacing;
  imSpacing[0] = 0.15;
  imSpacing[1] = 0.15;
  imSpacing[2] = 0.3;
  im->SetSpacing( imSpacing );
  im->Allocate();

  itk::ImageRegionIteratorWithIndex<ImageType> itIm( im,
    im->GetLargestPossibleRegion() );
  ImageType::PointType pnt;
  itIm.GoToBegin();
  double x, d, d2;
  while( !itIm.IsAtEnd() )
    {
    im->TransformIndexToPhysicalPoint( itIm.GetIndex(), pnt );
    x = pnt[0]+pnt[1]*0.2;
    if( itIm.GetIndex()[0] == itIm.GetIndex()[1]
      && itIm.GetIndex()[0] == 0 )
      {
      std::cout << "Slice index = "
        << ( itIm.GetIndex()[2]-index0[2] ) % 7 << std::endl;
      }
    switch( ( itIm.GetIndex()[2]-index0[2] ) % 7 )
      {
      default:
      case 0:
        {
        itIm.Set( spline.Value( x ) );
        break;
        }
      case 1:
        {
        itIm.Set( myFuncV->Value( x ) );
        break;
        }
      case 2:
        {
        itIm.Set( spline.ValueD( x ) );
        break;
        }
      case 3:
        {
        itIm.Set( myFuncD->Value( x ) );
        break;
        }
      case 4:
        {
        itIm.Set( spline.ValueD2( x ) );
        break;
        }
      case 5:
        {
        itIm.Set( spline.Curv( x ) );
        break;
        }
      case 6:
        {
        itIm.Set( spline.ValueJet( x, &d, &d2 ) );
        break;
        }
      }
    ++itIm;
    }

  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer imWriter = ImageWriterType::New();
  imWriter->SetFileName( argv[1] );
  imWriter->SetInput( im );
  imWriter->SetUseCompression( true );
  imWriter->Update();

  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer rndGen
    = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
  rndGen->Initialize( 1 );

  int failed = 0;
  for( unsigned int c=0; c<100; c++ )
    {
    x = rndGen->GetNormalVariate( 0.0, 0.5 );
    while( x > 1.4749 )
      {
      x = rndGen->GetNormalVariate( 0.0, 0.5 );
      }

    double xVal = 0;
    double x0 = x;
    if( !spline.Extreme( &x, &xVal ) )
      {
      std::cout << "Spline.Extreme() returned false." << std::endl;
      std::cout << "                x = " << x << std::endl;
      std::cout << "                xVal = " << xVal << std::endl;
      returnStatus = EXIT_FAILURE;
      ++failed;
      }
    else
      {
      bool err=false;
      if( vnl_math_abs( x - -1.4749 ) > 0.0001 )
        {
        std::cout << "Spline.Extreme() solution not ideal: x="
          << x << " != ideal=-1.4749" << std::endl;
        std::cout << "   Optimization started at x = " << x0 << std::endl;
        std::cout << "   Value at x0 = " << std::sin(x0) << std::endl;
        std::cout << "   Derivative at x0 = " << std::cos(x0) << std::endl;
        std::cout << "   Iteration = " << c << " of 100." << std::endl;
        for( double tx=x0-1; tx<x0+1; tx+=0.1 )
          {
          std::cout << "    tx0 = " << tx << " v = " << spline.Value( tx )
            << std::endl;
          }
        for( double tx=x-1; tx<x+1; tx+=0.1 )
          {
          std::cout << "    tx = " << tx << " v = " << spline.Value( tx )
            << std::endl;
          }
        returnStatus = EXIT_FAILURE;
        err = true;
        }
      if( vnl_math_abs( xVal - -0.827393 ) > 0.0001 )
        {
        std::cout << "Spline.Extreme() output not ideal: val="
          << xVal << " != ideal=1.0" << std::endl;
        std::cout << "   Optimization started at x = " << x0 << std::endl;
        std::cout << "   Iteration = " << c << " of 100." << std::endl;
        returnStatus = EXIT_FAILURE;
        err = true;
        }
      if( err )
        {
        ++failed;
        }
      }
    }

  delete myFunc;
  delete myFuncV;
  delete myFuncD;
  delete opt;

  std::cout << failed << " out of 100 optimizations failed." << std::endl;

  return returnStatus;
}
