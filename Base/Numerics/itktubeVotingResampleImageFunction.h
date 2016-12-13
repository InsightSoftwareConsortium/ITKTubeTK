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

#ifndef __itktubeVotingResampleImageFunction_h
#define __itktubeVotingResampleImageFunction_h

#include <itkInterpolateImageFunction.h>

namespace itk
{

namespace tube
{

/** \class VotingResampleImageFunction
 * \brief Linearly interpolate an image at specified positions.
 *
 * VotingResampleImageFunction linearly interpolates image intensity at
 * a non-integer pixel position. This class is templated
 * over the input image type and the coordinate representation type
 * ( e.g. float or double ).
 *
 * This function works for N-dimensional images.
 *
 * \warning This function work only for images with scalar pixel
 * types. For vector images use VectorVotingResampleImageFunction.
 *
 * \sa VectorVotingResampleImageFunction
 *
 * \ingroup ImageFunctions ImageInterpolators
 */
template< class TInputImage, class TCoordRep = float >
class VotingResampleImageFunction
  : public InterpolateImageFunction< TInputImage, TCoordRep >
{
public:
  /** Standard class typedefs. */
  typedef VotingResampleImageFunction                     Self;
  typedef InterpolateImageFunction<TInputImage,TCoordRep> Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( VotingResampleImageFunction, InterpolateImageFunction );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** OutputType typedef support. */
  typedef typename Superclass::OutputType OutputType;

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;

  /** RealType typedef support. */
  typedef typename Superclass::RealType RealType;

  /** Dimension underlying input image. */
  itkStaticConstMacro( ImageDimension, unsigned int, Superclass::ImageDimension );

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex( 
    const ContinuousIndexType & index ) const;

protected:
  VotingResampleImageFunction( void );
  ~VotingResampleImageFunction( void ) {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

private:
  VotingResampleImageFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long  m_Neighbors;

}; // End class VotingResampleImageFunction

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeVotingResampleImageFunction.hxx"
#endif

#endif // End !defined( __itktubeVotingResampleImageFunction_h )
