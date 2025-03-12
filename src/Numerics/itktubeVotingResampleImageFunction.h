/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef itktubeVotingResampleImageFunction_h
#define itktubeVotingResampleImageFunction_h

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
template <class TInputImage, class TCoordRep = float>
class VotingResampleImageFunction : public InterpolateImageFunction<TInputImage, TCoordRep>
{
public:
  /** Standard class type alias. */
  using Self = VotingResampleImageFunction;
  using Superclass = InterpolateImageFunction<TInputImage, TCoordRep>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro(VotingResampleImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** OutputType type alias support. */
  using OutputType = typename Superclass::OutputType;

  /** InputImageType type alias support. */
  using InputImageType = typename Superclass::InputImageType;

  /** SizeType type alias support. */
  using SizeType = typename Superclass::SizeType;

  /** RealType type alias support. */
  using RealType = typename Superclass::RealType;

  /** Dimension underlying input image. */
  static constexpr unsigned int ImageDimension = Superclass::ImageDimension;

  /** Index type alias support. */
  using IndexType = typename Superclass::IndexType;

  /** ContinuousIndex type alias support. */
  using ContinuousIndexType = typename Superclass::ContinuousIndexType;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType
  EvaluateAtContinuousIndex(const ContinuousIndexType & index) const override;

  SizeType
  GetRadius() const override
  {
    return m_Radius;
  }

protected:
  VotingResampleImageFunction(void);
  ~VotingResampleImageFunction(void) {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  VotingResampleImageFunction(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long m_Neighbors;
  SizeType                   m_Radius;

}; // End class VotingResampleImageFunction

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeVotingResampleImageFunction.hxx"
#endif

#endif // End !defined( __itktubeVotingResampleImageFunction_h )
