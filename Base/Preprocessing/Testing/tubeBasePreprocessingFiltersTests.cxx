/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBasicFiltersTests.cxx,v $
  Language:  C++
  Date:      $Date: 2009-11-02 15:26:53 $
  Version:   $Revision: 1.125 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <iostream>
#include "itkTestMain.h" 


void RegisterTests()
{
  REGISTER_TEST( tubeBasePreprocessingFiltersPrintTest );
  REGISTER_TEST( itkTubeEnhancingDiffusion2DImageFilterTest );
}

