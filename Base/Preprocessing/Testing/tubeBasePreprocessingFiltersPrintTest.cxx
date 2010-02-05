/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBasicFiltersPrintTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007-08-10 14:34:01 $
  Version:   $Revision: 1.20 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
#define ITK_LEAN_AND_MEAN

#include "itkVesselEnhancingDiffusion2DImageFilter.h"


int tubeBasePreprocessingFiltersPrintTest( int, char* [] )
{
  itk::VesselEnhancingDiffusion2DImageFilter< float, 2 >::Pointer 
       vesselEnahncingObj =
    itk::VesselEnhancingDiffusion2DImageFilter< float, 2 >::New();
  std::cout << "-------------VesselEnhancingDiffusion2DImageFilter" 
            << vesselEnahncingObj;

  return EXIT_SUCCESS;
}
