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

#include "itkGeometryObject.h"

using namespace itk;
using namespace itk::Geometry;

class GeometryObject3DTester : public GeometryObject<3>
{
public:
  typedef GeometryObject3DTester     Self;
  typedef GeometryObject<3>          Superclass; 
  
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer; 

  /** Method for creation through the object factory. */
  itkNewMacro( GeometryObject3DTester );
  itkTypeMacro( GeometryObject3DTester, GeometryObject );
};

class GeometryObject2DTester : public GeometryObject<2>
{
public:
  typedef GeometryObject2DTester     Self;
  typedef GeometryObject<2>          Superclass; 
  
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer; 

  /** Method for creation through the object factory. */
  itkNewMacro( GeometryObject2DTester );
  itkTypeMacro( GeometryObject2DTester, GeometryObject );
};

int itkGeometryObjectTest1(int itkNotUsed(argc), char* itkNotUsed(argc)[])
{
  GeometryObject3DTester::Pointer geom3d = GeometryObject3DTester::New();
  GeometryObject3DTester::Pointer geom2d = GeometryObject3DTester::New();
  
  return EXIT_SUCCESS;
}
