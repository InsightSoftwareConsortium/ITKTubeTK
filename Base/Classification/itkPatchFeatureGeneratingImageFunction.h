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
#ifndef __itkPatchFeatureGeneratingImageFunction_h
#define __itkPatchFeatureGeneratingImageFunction_h

#include <vector>

#include "itkImage.h"
#include "itkConstNeighborhoodIterator.h"

#include "itkFeatureGeneratingImageFunction.h"

namespace itk
{

/** \class PatchFeatureGeneratingImageFunction
 *
 */
template<class TInputImage, class TCoordRep = float>
class ITK_EXPORT PatchFeatureGeneratingImageFunction :
  public FeatureGeneratingImageFunction< TInputImage, TCoordRep >
{
public:

  /** Class typedefs **/
  typedef PatchFeatureGeneratingImageFunction       Self;
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
  typedef itk::ConstNeighborhoodIterator<InputImageType>             NeighborIterType;

  /** Run-time type information (and related methods). */
  itkTypeMacro( PatchFeatureGeneratingImageFunction, 
                FeatureGeneratingImageFunction );

  /** Patch New Macro. */
  itkNewMacro( Self );

  /** Constant for fetching the dimensions of the image. **/
  itkStaticConstMacro( ImageDimension, unsigned int,
                       Superclass::ImageDimension );
  
  void SetPriorImage( typename InputImageType::Pointer prior )
  {
    m_Prior = prior;
  }

  void SetWidth( size_t width )
  {
    m_PatchWidth = width;
  }

  /** Get the feature vector at an index for a given point **/
  virtual OutputType EvaluateAtIndex( const IndexType & index ) const;

protected:
  
  /** Default constructor */
  PatchFeatureGeneratingImageFunction();

  /** Default destructor */
  virtual ~PatchFeatureGeneratingImageFunction() {}

  typename InputImageType::Pointer m_Prior;

  size_t m_PatchWidth;

private:
  PatchFeatureGeneratingImageFunction( const Self& ); //pni
  void operator=( const Self& ); //purposely not implemented

}; // End class PatchFeatureGeneratingImageFunction

}// end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkPatchFeatureGeneratingImageFunction.txx"
#endif

#endif
