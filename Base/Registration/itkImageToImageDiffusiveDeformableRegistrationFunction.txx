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
#ifndef _itkImageToImageDiffusiveDeformableRegistrationFunction_txx
#define _itkImageToImageDiffusiveDeformableRegistrationFunction_txx

#include "itkImageToImageDiffusiveDeformableRegistrationFunction.h"

namespace itk
{ 
    
/**
 * Constructor
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::ImageToImageDiffusiveDeformableRegistrationFunction()
{
  m_TimeStep = 1.0;
  this->SetMovingImage(0);
  this->SetFixedImage(0);
}

/**
 * PrintSelf
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "TimeStep: " << m_TimeStep;
}

/**
  * Called at the beginning of each iteration
  */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::InitializeIteration()
{
}

/**
  * Called at the beginning of each iteration
  */
template < class TFixedImage, class TMovingImage, class TDeformationField >
typename ImageToImageDiffusiveDeformableRegistrationFunction
                                < TFixedImage, TMovingImage, TDeformationField >
::PixelType
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::ComputeUpdate(const NeighborhoodType &neighborhood,
                void *gd,
                const FloatOffsetType& itkNotUsed(offset))
{
  // Get the global data structure
  GlobalDataStruct * globalData = ( GlobalDataStruct * ) gd;
}

/**
  * Release the per-thread-global data
  */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
ImageToImageDiffusiveDeformableRegistrationFunction< TFixedImage,
                                                     TMovingImage,
                                                     TDeformationField >
::ReleaseGlobalDataPointer( void * gd ) const
{
  delete ( GlobalDataStruct * ) gd;
}






} // end namespace itk

#endif
