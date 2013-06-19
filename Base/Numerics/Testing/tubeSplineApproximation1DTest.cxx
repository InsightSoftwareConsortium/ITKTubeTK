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

#include "tubeBrentOptimizer1D.h"
#include "tubeSplineApproximation1D.h"
#include "tubeUserFunction.h"

#include <vcl_cmath.h>
#include <vnl/vnl_math.h>

#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

#include <cstdlib>
#include <iostream>

class MySA1DFunc : public tube::UserFunction< int, double >
{
private:
  double cVal;

public:
  MySA1DFunc( void )
    {
    cVal = 0;
    }
  const double & value( const int & x )
    {
    cVal = vcl_sin((double)x);
    std::cout << "s: x = " << x << " : v = " << cVal << std::endl;
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
  const double & value( const double & x )
    {
    cVal = vcl_sin((double)x);
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
  const double & value( const double & x )
    {
    cDeriv = vcl_cos((double)x);
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
  opt->smallDouble( epsilon );
  opt->searchForMin( true );
  opt->xStep( 0.01 );
  opt->tolerance( 0.0000001 );

  tube::SplineApproximation1D spline( myFunc, opt );

  int returnStatus = EXIT_SUCCESS;

  spline.clipEdge( true );

  spline.xMin( -3 );
  if( spline.xMin() != -3 )
    {
    std::cout << "xMin should be -3 and not " << spline.xMin() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  spline.xMax( 6 );
  if( spline.xMax() != 6 )
    {
    std::cout << "xMax should be 6 and not " << spline.xMax() << std::endl;
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
        << (itIm.GetIndex()[2]-index0[2]) % 7 << std::endl;
      }
    switch( (itIm.GetIndex()[2]-index0[2]) % 7 )
      {
      default:
      case 0:
        {
        itIm.Set( spline.value(x) );
        break;
        }
      case 1:
        {
        itIm.Set( myFuncV->value(x) );
        break;
        }
      case 2:
        {
        itIm.Set( spline.valueD(x) );
        break;
        }
      case 3:
        {
        itIm.Set( myFuncD->value(x) );
        break;
        }
      case 4:
        {
        itIm.Set( spline.valueD2(x) );
        break;
        }
      case 5:
        {
        itIm.Set( spline.curv(x) );
        break;
        }
      case 6:
        {
        itIm.Set( spline.valueJet(x, &d, &d2) );
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
  for(unsigned int c=0; c<100; c++)
    {
    x = rndGen->GetNormalVariate( 0.0, 0.5 );

    double xVal = 0;
    if( !spline.extreme( &x, &xVal ) )
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
          << x << "!= ideal=" << vnl_math::pi/2 << std::endl;
        returnStatus = EXIT_FAILURE;
        err = true;
        }
      if( vnl_math_abs( xVal - -0.827393 ) > 0.0001 )
        {
        std::cout << "Spline.Extreme() output not ideal: val="
          << xVal << " != ideal=1.0" << std::endl;
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
