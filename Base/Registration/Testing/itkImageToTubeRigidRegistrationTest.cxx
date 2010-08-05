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


//#include "itkImageToTubeRigidRegistration.h"
#include "itkDefaultDynamicMeshTraits.h"
#include "itkSpatialObjectReader.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int itkImageToTubeRigidRegistrationTest(int argc, char* argv [] )
{
  if ( argc < 2 )
    {
    std::cerr << "Missing Parameters: " 
              << argv[0]
              << " Input_Vessel" << std::endl;
//              << " Edge_Enhanced_Output_Image [Sigma] [Alpha] [ContrastParameter]"
//              << " [TimeStep] [NumberOfIterations]" << std::endl; 
    return EXIT_FAILURE;
    }
 
 
  typedef itk::DefaultDynamicMeshTraits< float , 3, 3>    MeshTraitType;
  typedef itk::SpatialObjectReader<3,float,MeshTraitType> SpatialObjectReaderType;

  SpatialObjectReaderType::Pointer vesselReader = SpatialObjectReaderType::New();
  vesselReader->SetFileName(argv[1]);
  
  try
    {
    vesselReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;

}

