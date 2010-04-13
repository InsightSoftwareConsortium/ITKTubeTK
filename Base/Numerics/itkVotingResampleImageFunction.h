/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVotingResampleImageFunction.h,v $
  Language:  C++
  Date:      $Date: 2007/01/30 20:56:09 $
  Version:   $Revision: 1.31 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVotingResampleImageFunction_h
#define __itkVotingResampleImageFunction_h

#include "itkInterpolateImageFunction.h"

namespace itk
{

/** \class VotingResampleImageFunction
 * \brief Linearly interpolate an image at specified positions.
 *
 * VotingResampleImageFunction linearly interpolates image intensity at
 * a non-integer pixel position. This class is templated
 * over the input image type and the coordinate representation type 
 * (e.g. float or double).
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
class ITK_EXPORT VotingResampleImageFunction : 
  public InterpolateImageFunction<TInputImage,TCoordRep> 
{
public:
  /** Standard class typedefs. */
  typedef VotingResampleImageFunction                     Self;
  typedef InterpolateImageFunction<TInputImage,TCoordRep> Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(VotingResampleImageFunction, InterpolateImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** OutputType typedef support. */
  typedef typename Superclass::OutputType OutputType;

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;

  /** RealType typedef support. */
  typedef typename Superclass::RealType RealType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

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
  VotingResampleImageFunction();
  ~VotingResampleImageFunction(){};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  VotingResampleImageFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long  m_Neighbors;  

};

} // end namespace itk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_VotingResampleImageFunction(_, EXPORT, x, y) \
  namespace itk { \
  _(2(class EXPORT VotingResampleImageFunction< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef VotingResampleImageFunction< ITK_TEMPLATE_2 x > \
                                              VotingResampleImageFunction##y; } \
  }

#if ITK_TEMPLATE_TXX
# include "itkVotingResampleImageFunction.txx"
#endif

#endif
