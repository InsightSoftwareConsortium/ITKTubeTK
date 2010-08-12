/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAnisotropicCoherenceEnhancingDiffusionImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007/06/20 16:03:23 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include "itkTubeToTubeTransformFilter.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkEuler3DTransform.h"
#include "itkMath.h"

int itkTubeToTubeTransformFilterTest(int argc, char* argv [] )
{

  if ( argc < 3 )
    {
    std::cerr << "Missing Parameters: " 
              << argv[0]
              << " Input_Vessel " << "Input_Image " << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::GroupSpatialObject<3>                      TubeNetType;
  typedef itk::SpatialObjectReader<3>                     TubeNetReaderType;
  typedef itk::SpatialObjectWriter<3>                     TubeNetWriterType;  
  typedef itk::Euler3DTransform<double>                   TransformType; 
  typedef itk::TubeToTubeTransformFilter<TransformType,3> TubeTransformFilterType;

  // read in vessel
  TubeNetReaderType::Pointer reader = TubeNetReaderType::New();
  reader->SetFileName(argv[1]);
  
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // generate transform
  TransformType::Pointer transform = TransformType::New();

  itk::Vector<double,3> CoR;
  CoR[0] = 0.0;
  CoR[1] = 0.0;
  CoR[2] = 0.0;

  itk::Vector<double,3> rotation;
  rotation[0] = 0;//-0.5/itk::Math::one_over_pi;
  rotation[1] = 0;//-0.5/itk::Math::one_over_pi;
  rotation[2] = 0.0;

  itk::Vector<double,3> translation;
  translation[0] = 0.0;
  translation[1] = 0.0;
  translation[2] = 0.0;

  double ca=cos(rotation[0]);
  double sa=sin(rotation[0]);
  double cb=cos(rotation[1]);
  double sb=sin(rotation[1]);
  double cg=cos(rotation[2]);
  double sg=sin(rotation[2]);

  itk::Matrix<double,3,3> rotationMatrix;
  rotationMatrix[0][0] = ca*cb;
  rotationMatrix[0][1] = ca*sb*sg-sa*cg;
  rotationMatrix[0][2] = ca*sb*cg+sa*sg;
  rotationMatrix[1][0] = sa*cb;
  rotationMatrix[1][1] = sa*sb*sg+ca*cg;
  rotationMatrix[1][2] = sa*sb*cg-ca*sg;
  rotationMatrix[2][0] = -sb;
  rotationMatrix[2][1] = cb*sg;
  rotationMatrix[2][2] = cb*cg;

  transform->SetRotationMatrix(rotationMatrix);

  //transform->SetRotation(rotation[0], rotation[1], rotation[2]);
  
  itk::Vector<double,3> rotOffset = -(transform->GetRotationMatrix()*CoR);

  for(unsigned int i=0;i<3;i++)
    {
    rotOffset[i] += translation[i]+CoR[i];
    }
  
  transform->SetOffset(rotOffset);

  // create transform filter
  TubeTransformFilterType::Pointer transformFilter = TubeTransformFilterType::New();
  transformFilter->SetInput(reader->GetGroup());
  transformFilter->SetScale(1.0);
  transformFilter->SetTransform(transform);
  
  try
    {
    transformFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }
  
  // write vessel
  TubeNetWriterType::Pointer writer = TubeNetWriterType::New();
  writer->SetFileName(argv[2]);
  writer->SetInput(transformFilter->GetOutput());

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

