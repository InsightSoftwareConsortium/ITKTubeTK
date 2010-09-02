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
#ifndef __itkStandardFeatureGeneratingImageFunction_txx
#define __itkStandardFeatureGeneratingImageFunction_txx

#include "itkStandardFeatureGeneratingImageFunction.h"

#include "itkRidgeExtractor.h"
#include "itkJointHistogramImageFunction.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkDivideImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "vnl/vnl_math.h"
#include "math.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>

namespace itk
{

/**
 * Set the input Image
 */
template <class TInputImage, class TCoordRep>
StandardFeatureGeneratingImageFunction<TInputImage,TCoordRep>
::StandardFeatureGeneratingImageFunction()
{
  this->m_Features.push_back("RidgenessDiffSmall");
  this->m_Features.push_back("RidgenessDiffMedium");
  this->m_Features.push_back("RidgenessDiffLarge");


}

template <class TInputImage, class TCoordRep>
typename StandardFeatureGeneratingImageFunction<TInputImage,TCoordRep>::OutputType
StandardFeatureGeneratingImageFunction<TInputImage,TCoordRep>
::EvaluateAtIndex( const IndexType & index ) const
{
  OutputType x( this->m_Features.size() );
  PointType pnt;
  this->GetInputImage()->TransformIndexToPhysicalPoint( index,
                                                        pnt );
  x[0] = pnt[0];
  x[1] = pnt[1];
  x[2] = this->GetInputImage()->GetPixel( index );
  return x;
}

}

#endif
