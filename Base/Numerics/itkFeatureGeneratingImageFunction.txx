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
#ifndef __itkFeatureGeneratingImageFunction_txx
#define __itkFeatureGeneratingImageFunction_txx

#include "itkFeatureGeneratingImageFunction.h"

namespace itk
{

/**
 * Set the input Image
 */
template <class TInputImage, class TCoordRep>
FeatureGeneratingImageFunction<TInputImage,TCoordRep>
::FeatureGeneratingImageFunction()
{
}

template <class TInputImage, class TCoordRep>
void
FeatureGeneratingImageFunction<TInputImage,TCoordRep>
::SetInputImage( const InputImageType * ptr )
{
  this->Superclass::SetInputImage( ptr );
}

template <class TInputImage, class TCoordRep>
typename FeatureGeneratingImageFunction<TInputImage,TCoordRep>::OutputType
FeatureGeneratingImageFunction<TInputImage,TCoordRep>
::EvaluateAtIndex( const IndexType & index ) const
{
  vnl_vector<double> x(3);
  PointType pnt;
  this->GetInputImage()->TransformIndexToPhysicalPoint( index, 
                                                        pnt );
  x[0] = pnt[0];
  x[1] = pnt[1];
  x[2] = this->GetInputImage()->GetPixel( index );
  return x;
}

template <class TInputImage, class TCoordRep>
void
FeatureGeneratingImageFunction<TInputImage,TCoordRep>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );
}

} // end namespace itk

#endif
