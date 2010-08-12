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

#include "../itkSplineND.h"
#include "../itkOptBrent1D.h"
#include "../itkSplineApproximation1D.h"
#include "../itkUserFunc.h"

class MySANDFunc:
  public itk::UserFunc< vnl_vector<int>, double > 
  {
  private:
    double cVal;
  public:
    MySANDFunc( )
      {
      cVal = 0;
      };
    const double & value( const vnl_vector<int> & x )
      {
      cVal = vcl_sin((double)x[0]);
      cVal += vcl_cos((double)x[1]);
      std::cout << "s: x = " << x[0] << ", " << x[1] 
        << " : v = " << cVal << std::endl;
      return cVal;
      };
  };

int itkSplineNDTest( int argc, char *argv[] )
{
  if( argc != 2 )
    {
    std::cout << "usage: run <outImFile>" << std::endl;
    return EXIT_FAILURE;
    }

  double epsilon = 0.000001;

  MySANDFunc * myFunc = new MySANDFunc();

  itk::SplineApproximation1D * spline1D = new itk::SplineApproximation1D();

  itk::OptBrent1D * opt = new itk::OptBrent1D( );
  opt->smallDouble( epsilon );

  itk::SplineND spline( 2, myFunc, spline1D, opt );

  int returnStatus = EXIT_SUCCESS;

  vnl_vector<int> xMin(2, -3);
  spline.xMin( xMin );
  if( spline.xMin()[0] != -3 )
    {
    std::cout << "xMin should be -3 and not " << spline.xMin() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  vnl_vector<int> xMax(2, 3);
  spline.xMax( xMax );
  if( spline.xMax()[0] != 3 )
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
  vnl_vector<double> x(2);
  vnl_vector<double> d(2);
  vnl_vector<int> di(2);
  di.fill( 0 );
  di[1] = 1;
  vnl_matrix<double> m(2,2);
  vnl_vector<double> d2(2);
  while( !itIm.IsAtEnd() )
    {
    im->TransformIndexToPhysicalPoint( itIm.GetIndex(), pnt );
    x[0] = pnt[0];
    x[1] = pnt[1];
    switch( vnl_math_abs(itIm.GetIndex()[2]) % 5 )
      {
      default:
      case 0:
        {
        itIm.Set( spline.value(x) );
        break;
        }
      case 1:
        {
        itIm.Set( spline.valueD(x, di) );
        break;
        }
      case 2:
        {
        d = spline.valueD(x);
        itIm.Set( d[0]*d[0] + d[1]*d[1] );
        break;
        }
      case 3:
        {
        m = spline.hessian(x);
        itIm.Set( d[0]*d[0] + d[1]*d[1] );
        break;
        }
      case 4:
        {
        itIm.Set( spline.valueJet(x, d, m) );
        break;
        }
      case 5:
        {
        itIm.Set( spline.valueVDD2(x, d, d2) );
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

  return returnStatus;
}

