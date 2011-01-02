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

#include <sstream>

#include "itkFeatureGeneratingImageFunction.h"

namespace itk
{

/**
 * Set the input Image
 */
template <class TInputImage, class TCoordRep>
FeatureGeneratingImageFunction<TInputImage,TCoordRep>
::FeatureGeneratingImageFunction() :
  m_Features()
{
  m_Features.push_back("x");
  m_Features.push_back("y");
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
  OutputType x( m_Features.size() );
  PointType pnt;
  this->GetInputImage()->TransformIndexToPhysicalPoint( index,
                                                        pnt );
  x[0] = pnt[0];
  x[1] = pnt[1];
  return x;
}

template <class TInputImage, class TCoordRep>
std::string
FeatureGeneratingImageFunction<TInputImage,TCoordRep>
::EvaluateToStringAtIndex( const IndexType & index ) const
{
  OutputType out = this->EvaluateAtIndex( index );
  std::stringstream ss;
  typename OutputType::const_iterator itr;
  typename OutputType::const_iterator itrPlusOne;
  for( itr = out.begin(), itrPlusOne = out.begin(); itr != out.end(); ++itr )
    {
    ++itrPlusOne;
    if( itrPlusOne != out.end() )
      {
      ss << *itr << ",";
      }
    else
      {
      ss << *itr;
      }
    }
  return ss.str();
}

template <class TInputImage, class TCoordRep>
const typename FeatureGeneratingImageFunction<TInputImage,TCoordRep>::FeatureListType&
FeatureGeneratingImageFunction<TInputImage,TCoordRep>
::GetFeatureLabels() const
{
  return m_Features;
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
