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
#ifndef __itkStandardFeatureGeneratingImageFunction_h
#define __itkStandardFeatureGeneratingImageFunction_h

#include <vector>

#include "itkImageFunction.h"
#include "itkOrientedImage.h"
#include "itkFeatureGeneratingImageFunction.h"

namespace itk
{

/** \class StandardFeatureGeneratingImageFunction
 *
 */
template<class TInputImage, class TCoordRep = float>
class ITK_EXPORT StandardFeatureGeneratingImageFunction :
  public FeatureGeneratingImageFunction< TInputImage, TCoordRep >
{
public:

  /** Class typedefs **/
  typedef StandardFeatureGeneratingImageFunction       Self;
  typedef FeatureGeneratingImageFunction<TInputImage,TCoordRep>  
                                                       Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;
  typedef typename Superclass::InputImageType          InputImageType;
  typedef typename TInputImage::PixelType              PixelType;
  typedef typename Superclass::PointType               PointType;
  typedef typename Superclass::IndexType               IndexType;
  typedef typename Superclass::ContinuousIndexType     ContinuousIndexType;
  typedef typename Superclass::OutputType              OutputType;
  typedef typename Superclass::FeatureListType         FeatureListType;

  /** Run-time type information (and related methods). */
  itkTypeMacro( StandardFeatureGeneratingImageFunction, 
                FeatureGeneratingImageFunction );

  /** Standard New Macro. */
  itkNewMacro( Self );

  /** Constant for fetching the dimensions of the image. **/
  itkStaticConstMacro( ImageDimension, unsigned int,
                       Superclass::ImageDimension );

  /** Get the feature vector at an index for a given point **/
  virtual OutputType EvaluateAtIndex( const IndexType & index ) const;

protected:
  
  /** Default constructor */
  StandardFeatureGeneratingImageFunction();

  /** Default destructor */
  ~StandardFeatureGeneratingImageFunction() {}

private:
  StandardFeatureGeneratingImageFunction( const Self& ); //pni
  void operator=( const Self& ); //purposely not implemented

}; // End class StandardFeatureGeneratingImageFunction

}// end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkStandardFeatureGeneratingImageFunction.txx"
#endif

#endif
