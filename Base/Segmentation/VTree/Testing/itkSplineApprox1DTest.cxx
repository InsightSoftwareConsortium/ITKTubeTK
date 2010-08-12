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
#include <cstdlib>
#include <iostream>

#include <vnl/vnl_math.h>
#include <vcl_cmath.h>

#include "itkMacro.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"

#include "../itkOptBrent1D.h"
#include "../itkSplineApproximation1D.h"
#include "../itkUserFunc.h"

class MySA1DFunc:
  public itk::UserFunc< int, double > 
  {
  private:
    double cVal;
  public:
    MySA1DFunc( )
      {
      cVal = 0;
      };
    const double & value( const int & x )
      {
      cVal = vcl_sin((double)x);
      std::cout << "s: x = " << x << " : v = " << cVal << std::endl;
      return cVal;
      };
  };

int itkSplineApprox1DTest( int argc, char *argv[] )
{
  if( argc != 2 )
    {
    std::cout << "usage: run <outImFile>" << std::endl;
    return EXIT_FAILURE;
    }

  double epsilon = 0.000001;

  MySA1DFunc * myFunc = new MySA1DFunc();

  itk::OptBrent1D * opt = new itk::OptBrent1D( );
  opt->smallDouble( epsilon );
  opt->searchForMin( true );
  opt->xStep( 0.01 );
  opt->tolerance( 0.0000001 );

  itk::SplineApproximation1D spline( myFunc, opt );

  int returnStatus = EXIT_SUCCESS;

  spline.xMin( -3 );
  if( spline.xMin() != -3 )
    {
    std::cout << "xMin should be -3 and not " << spline.xMin() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  spline.xMax( 3 );
  if( spline.xMax() != 3 )
    {
    std::cout << "xMax should be 3 and not " << spline.xMax() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  typedef itk::Image< float, 3 >  ImageType;

  ImageType::Pointer im = ImageType::New( );
  ImageType::RegionType imRegion;
  ImageType::SizeType imSize;
  imSize[0] = 40;
  imSize[1] = 40;
  imSize[2] = 20;
  imRegion.SetSize( imSize );
  ImageType::IndexType index0;
  index0[0] = -20;
  index0[1] = -20;
  index0[2] = -10;
  imRegion.SetIndex( index0 );
  im->SetRegions( imRegion );
  ImageType::SpacingType imSpacing;
  imSpacing[0] = 0.1;
  imSpacing[1] = 0.1;
  imSpacing[2] = 0.2;
  im->SetSpacing( imSpacing );
  im->Allocate( );

  itk::ImageRegionIteratorWithIndex<ImageType> itIm( im,
    im->GetLargestPossibleRegion( ) );
  ImageType::PointType pnt;
  itIm.GoToBegin();
  double x, d, d2;
  while( !itIm.IsAtEnd() )
    {
    im->TransformIndexToPhysicalPoint( itIm.GetIndex(), pnt );
    x = pnt[0];
    switch( vnl_math_abs(itIm.GetIndex()[1]) % 5 )
      {
      default:
      case 0:
        {
        itIm.Set( spline.value(x) );
        break;
        }
      case 1:
        {
        itIm.Set( spline.valueD(x) );
        break;
        }
      case 2:
        {
        itIm.Set( spline.valueD2(x) );
        break;
        }
      case 3:
        {
        itIm.Set( spline.curv(x) );
        break;
        }
      case 4:
        {
        itIm.Set( spline.valueJet(x, &d, &d2) );
        break;
        }
      };
    ++itIm;
    }

  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer imWriter = ImageWriterType::New( );
  imWriter->SetFileName( argv[1] );
  imWriter->SetInput( im );
  imWriter->Update( );

  x = 0;
  double xVal;
  if( !spline.extreme( &x, &xVal ) )
    {
    std::cout << "Spline.Extreme() returned false." << std::endl;
    std::cout << "                x = " << x << std::endl;
    std::cout << "                xVal = " << xVal << std::endl;
    returnStatus = EXIT_FAILURE;
    }
  else
    {
    if( vnl_math_abs( x - -1.4749 ) > 0.0001 )
      {
      std::cout << "Spline.Extreme() solution not ideal: x=" 
        << x << "!= ideal=" << vnl_math::pi/2 << std::endl;
      returnStatus = EXIT_FAILURE;
      }
    if( vnl_math_abs( xVal - -0.827393 ) > 0.0001 )
      {
      std::cout << "Spline.Extreme() output not ideal: val=" 
        << xVal << " != ideal=1.0" << std::endl;
      returnStatus = EXIT_FAILURE;
      }
    }



  return returnStatus;
}

