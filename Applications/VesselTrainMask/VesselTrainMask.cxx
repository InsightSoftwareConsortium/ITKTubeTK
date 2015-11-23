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
#include <iostream>
#include <sstream>

#include "itkTimeProbesCollectorBase.h"
#include "tubeMessage.h"

#include "tubeMacro.h"
#include "metaScene.h"

#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkGroupSpatialObject.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "VesselTrainMaskCLP.h"

#include <vtkCubeSource.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>

// Description:
// Get the ComponentType and dimension of the image
void GetImageInformation(std::string fileName,
  itk::ImageIOBase::IOComponentType &componentType,
  unsigned int & dimension)
{
  itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO(fileName.c_str(),
    itk::ImageIOFactory::ReadMode);
  if (!imageIO)
    {
    std::cerr << "NO IMAGEIO WAS FOUND" << std::endl;
    dimension = 0;
    return;
  }

  // Now that we found the appropriate ImageIO class, ask it to
  // read the meta data from the image file.
  imageIO->SetFileName(fileName.c_str());
  imageIO->ReadImageInformation();

  componentType = imageIO->GetComponentType();
  dimension = imageIO->GetNumberOfDimensions();
}

// Main
int main( int argc, char * argv[] )
{
  try
    {
    PARSE_ARGS;
    }
  catch( const std::exception & err )
    {
    tube::ErrorMessage( err.what() );
    return EXIT_FAILURE;
    }
  PARSE_ARGS;

}
