/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeTubePointWeightsCalculator_h
#define __itktubeTubePointWeightsCalculator_h

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <itkSpatialObject.h>

namespace itk
{

namespace tube
{

/**
 *  This class computes scalar weights for every point in a tube tree
 *  based on
 *  the radius at that point.
 *
 *  \tparam TTubeTreeSpatialObject input tube tree spatial object type.
 *  \tparam TPointWeightFunction type of the function used to compute the
 *  weights.
 *  \tparam TResolutionsWeights type of the output scalar resolution
 *  weights.
 */
template< unsigned int VDimension, class TTubeSpatialObject,
  class TPointWeightFunction,
  class TPointWeights >
class TubePointWeightsCalculator : public Object
{
public:
  /** Standard class typedefs. */
  typedef TubePointWeightsCalculator< VDimension,
                                           TTubeSpatialObject,
                                           TPointWeightFunction,
                                           TPointWeights > Self;
  typedef Object                                           Superclass;
  typedef SmartPointer< Self >                             Pointer;
  typedef SmartPointer< const Self >                       ConstPointer;

  static const unsigned int Dimension = VDimension;

  typedef SpatialObject< Dimension > TubeTreeSpatialObjectType;
  typedef TTubeSpatialObject         TubeSpatialObjectType;
  typedef TPointWeightFunction       PointWeightFunctionType;
  typedef TPointWeights              PointWeightsType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( TubePointWeightsCalculator, Object );

  /** Compute the resolutions weights on the input TubeTreeSpatialObject
   * using the
   * PointWeightFunction. */
  void Compute( void );

  /** Set the input TubeTreeSpatialObject. */
  itkSetObjectMacro( TubeTreeSpatialObject, TubeTreeSpatialObjectType );
  itkGetConstObjectMacro( TubeTreeSpatialObject, TubeTreeSpatialObjectType );

  /** Set/Get the function used to determine the resolution weights.
   * This function
   *  takes a tube point as an input and outputs a weight for that point. */
  itkSetObjectMacro( PointWeightFunction, PointWeightFunctionType );
  itkGetConstObjectMacro( PointWeightFunction, PointWeightFunctionType );

  /** Get the output resolution weights. */
  itkGetConstReferenceMacro( PointWeights, PointWeightsType );

protected:
  TubePointWeightsCalculator( void );
  virtual ~TubePointWeightsCalculator( void ) {}

  void PrintSelf( std::ostream& os, Indent indent ) const;

  typename PointWeightFunctionType::ConstPointer m_PointWeightFunction;

  PointWeightsType m_PointWeights;

private:
  TubePointWeightsCalculator( const Self& ); //purposely not implemented
  void operator=( const Self& );            //purposely not implemented

  typename TubeTreeSpatialObjectType::ConstPointer m_TubeTreeSpatialObject;

}; // End class TubePointWeightsCalculator

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTubePointWeightsCalculator.hxx"
#endif

#endif // End !defined( __itktubeTubePointWeightsCalculator_h )
