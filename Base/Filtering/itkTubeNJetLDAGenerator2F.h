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
#ifndef __itkTubeNJetLDAGenerator2F_h
#define __itkTubeNJetLDAGenerator2F_h

#include <vector>

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include "itkImage.h"

#include "itkTubeLDAGenerator.h"
#include "itkTubeNJetLDAGenerator.h"

namespace itk
{

namespace tube
{

class NJetLDAGenerator2F
: public NJetLDAGenerator< Image< float, 2 >,
  Image< unsigned char, 2 > >
{
public:

  typedef NJetLDAGenerator2F                   Self;
  typedef NJetLDAGenerator< Image< float, 2 >,
    Image< unsigned char, 2 > >        Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkTypeMacro( NJetLDAGenerator2F, NJetLDAGenerator );

  itkNewMacro( Self );

  //
  // Custom Typedefs
  //
  typedef Image< float, 2 >                     ImageT;
  typedef ImageT                                        ImageType;
  typedef std::vector< ImageType::Pointer >             ImageListType;

  itkStaticConstMacro( ImageDimension, unsigned int,
    ImageT::ImageDimension );

  typedef Image< unsigned char, 2 >    LabelmapT;
  typedef LabelmapT                            MaskImageType;

  typedef double                               FeatureType;
  typedef vnl_vector< FeatureType >            FeatureVectorType;

  typedef int                                  ObjectIdType;
  typedef std::vector< ObjectIdType >          ObjectIdListType;

  typedef vnl_vector< double >                 ObjectMeanType;
  typedef std::vector< ObjectMeanType >        ObjectMeanListType;

  typedef vnl_matrix< double >                 ObjectCovarianceType;
  typedef std::vector< ObjectCovarianceType >  ObjectCovarianceListType;

  typedef vnl_vector< double >                 LDAValuesType;
  typedef vnl_vector< double >                 LDAVectorType;
  typedef vnl_matrix< double >                 LDAMatrixType;

  typedef std::vector< double >                NJetScalesType;

  typedef itk::Image< float, ImageDimension >  LDAImageType;
  typedef std::vector< LDAImageType::Pointer > LDAImageListType;

protected:

  NJetLDAGenerator2F( void );
  virtual ~NJetLDAGenerator2F( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  NJetLDAGenerator2F( const Self & );    // Purposely not implemented
  void operator = ( const Self & );      // Purposely not implemented

};

}

}

#endif
