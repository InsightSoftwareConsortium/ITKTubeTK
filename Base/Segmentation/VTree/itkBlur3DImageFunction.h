/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkBlur3DImageFunction.h,v $
  Language:  C++
  Date:      $Date: 2005/10/05 16:47:12 $
  Version:   $Revision: 1.6 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkBlur3DImageFunction_h
#define _itkBlur3DImageFunction_h

#include <itkIndex.h>
#include <itkImageFunction.h>

namespace itk
{

/**
 * \class Blur3DImageFunction
 * \brief Calculate the gaussian blurred value at point 
 *        given a scale and extent of the gaussian.
 * This class is templated over the input image type.
 *
 */
template <class TInputImage>
class ITK_EXPORT Blur3DImageFunction :
  public ImageFunction< TInputImage, double, double >
{
public:
  /**
   * Standard "Self" typedef */
  typedef Blur3DImageFunction                          Self;
  typedef ImageFunction<TInputImage, double, double>   Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  itkTypeMacro(Blur3DImageFunction, ImageFunction);

  itkNewMacro(Self);

  /**
   * InputImageType typedef support. */
  typedef TInputImage                                InputImageType;
  typedef typename InputImageType::SpacingType       SpacingType;

  /**
   * IndexType typedef support. */
  typedef typename InputImageType::IndexType         IndexType;
  typedef typename Superclass::ContinuousIndexType   ContinuousIndexType;

  /**
   * Dimension of the underlying image. */
  itkStaticConstMacro(ImageDimension, unsigned int, 3 ); 

  /**
   * Point typedef support. */
  typedef typename Superclass::PointType             PointType;

  /**
   * Set the input image. */
  virtual void SetInputImage( const InputImageType * ptr ); 

  /**
   * Evalulate the function at specified point */
  virtual double Evaluate( const PointType& point ) const;

  /** Evaluate the function at specified Index position. */
  virtual double EvaluateAtIndex( const IndexType & index ) const;

  /** Evaluate the function at specified ContinousIndex position. */
  virtual double EvaluateAtContinuousIndex( 
    const ContinuousIndexType & index ) const ;

  /**
   * */
  double EvaluateLinearInterpolate(const ContinuousIndexType & point) const;

  /**
   * Set the Scale */
  void SetScale(double scale);

  /**
   * Get the Scale */
   itkGetMacro( Scale, double);

  /**
   * Set the Extent */
  void SetExtent(double extent);

  /**
   * Get the Extent */
  itkGetMacro( Extent, double);

  /**
   * Get the Spacing */
  itkGetMacro( Spacing, SpacingType );
 
  /**
   * Interpret the sigma value to be in terms of x-spacing */
  void SetUseRelativeSpacing( bool useRelativeSpacing);
 
  /**
   * Get the Spacing */
  itkGetMacro( UseRelativeSpacing, bool );
 
protected:

  Blur3DImageFunction();
  virtual ~Blur3DImageFunction(){};

  void PrintSelf(std::ostream& os, Indent indent) const;

  void RecomputeKernel( void );

private:

  Blur3DImageFunction( const Self& );
  void operator=( const Self& );

  typedef std::list<double>                       KernelWeightsListType;
  typedef std::list< Index<3> >                   KernelXListType;

  bool                    m_Debug;

  Size<itkGetStaticConstMacro(ImageDimension)>    m_ImageSize;
  bool                    m_UseRelativeSpacing;
  SpacingType             m_Spacing;
  SpacingType             m_OriginalSpacing;
  double                  m_Scale;
  double                  m_Extent;
  KernelWeightsListType   m_KernelWeights;
  KernelXListType         m_KernelX;
  IndexType               m_KernelMin;
  IndexType               m_KernelMax;
  double                  m_KernelTotal;
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlur3DImageFunction.txx"
#endif

#endif

