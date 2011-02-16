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
#ifndef __itkAnisotropicDiffusiveSparseRegistrationFilter_txx
#define __itkAnisotropicDiffusiveSparseRegistrationFilter_txx

#include "itkAnisotropicDiffusiveSparseRegistrationFilter.h"

namespace itk
{

/**
 * Constructor
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
AnisotropicDiffusiveSparseRegistrationFilter
< TFixedImage, TMovingImage, TDeformationField >
::AnisotropicDiffusiveSparseRegistrationFilter()
{
}

/**
 * PrintSelf
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

/**
 * Setup the pointers for the deformation component images
 */
template < class TFixedImage, class TMovingImage, class TDeformationField >
void
AnisotropicDiffusiveSparseRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::InitializeDeformationComponentAndDerivativeImages()
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );

  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  OutputImagePointer output = this->GetOutput();

  // Setup pointers to the deformation component images - we have the
  // TANGENTIAL component, which is the entire deformation field, and the
  // NORMAL component, which is the deformation vectors projected onto their
  // normals
  // There are four terms that share these two deformation component images
  this->SetDeformationComponentImage( SMOOTH_TANGENTIAL, this->GetOutput() );
  this->SetDeformationComponentImage( PROP_TANGENTIAL, this->GetOutput() );

  DeformationFieldPointer normalDeformationField = DeformationFieldType::New();
  this->AllocateSpaceForImage( normalDeformationField, output );
  this->SetDeformationComponentImage( SMOOTH_NORMAL, normalDeformationField );
  this->SetDeformationComponentImage( PROP_NORMAL, normalDeformationField );

  // Setup the first and second order deformation component images
  // The two TANGENTIAL and two NORMAL components share images, so we allocate
  // images only for the SMOOTH terms and set the pointers for the PROP terms
  // to zero.  They will be pointed to the correct images when the images are
  // actually filled in.
  for( int i = 0; i < this->GetNumberOfTerms(); i++ )
    {
    ScalarDerivativeImageArrayType firstOrderArray;
    TensorDerivativeImageArrayType secondOrderArray;
    if( i % 2 == 0 )
      {
      for( int j = 0; j < ImageDimension; j++ )
        {
        firstOrderArray[j] = ScalarDerivativeImageType::New();
        this->AllocateSpaceForImage( firstOrderArray[j], output );
        secondOrderArray[j] = TensorDerivativeImageType::New();
        this->AllocateSpaceForImage( secondOrderArray[j], output );
        }
      }
    else
      {
      for( int j = 0; j < ImageDimension; j++ )
        {
        firstOrderArray[j] = 0;
        secondOrderArray[j] = 0;
        }
      }
    this->SetDeformationComponentFirstOrderDerivativeArray( i,
                                                            firstOrderArray );
    this->SetDeformationComponentSecondOrderDerivativeArray( i,
                                                             secondOrderArray );
    }
}






} // end namespace itk

#endif
