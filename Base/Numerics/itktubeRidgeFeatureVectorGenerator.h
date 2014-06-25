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

#ifndef __itktubeRidgeFeatureVectorGenerator_h
#define __itktubeRidgeFeatureVectorGenerator_h

#include "itktubeFeatureVectorGenerator.h"

#include <itkImage.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <vector>

namespace itk
{

namespace tube
{

template< class TImage >
class RidgeFeatureVectorGenerator
  : public FeatureVectorGenerator< Image< typename TImage::PixelType,
                                          TImage::ImageDimension > >
{
public:

  typedef RidgeFeatureVectorGenerator          Self;
  typedef FeatureVectorGenerator< TImage >     Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkTypeMacro( RidgeFeatureVectorGenerator, FeatureVectorGenerator );

  itkNewMacro( Self );

  //
  itkStaticConstMacro( ImageDimension, unsigned int,
    TImage::ImageDimension );

  typedef typename Superclass::FeatureValueType      FeatureValueType;

  typedef typename Superclass::FeatureVectorType     FeatureVectorType;

  typedef typename Superclass::ImageType             ImageType;

  typedef typename Superclass::IndexType             IndexType;

  typedef std::vector< double >                      RidgeScalesType;

  //
  virtual unsigned int GetNumberOfFeatures( void ) const;

  void SetScales( const RidgeScalesType & scales );
  const RidgeScalesType & GetScales( void ) const;

  virtual FeatureVectorType GetFeatureVector(
    const IndexType & indx ) const;

  virtual FeatureValueType GetFeatureVectorValue( const IndexType & indx,
    unsigned int fNum ) const;

protected:

  RidgeFeatureVectorGenerator( void );
  virtual ~RidgeFeatureVectorGenerator( void );

  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  // Purposely not implemented
  RidgeFeatureVectorGenerator( const Self & );
  void operator = ( const Self & );      // Purposely not implemented

  RidgeScalesType                    m_Scales;

}; // End class RidgeFeatureVectorGenerator

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeRidgeFeatureVectorGenerator.hxx"
#endif

#endif // End !defined(__itktubeRidgeFeatureVectorGenerator_h)
