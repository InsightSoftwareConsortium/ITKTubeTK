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
#ifndef __itkBlur3DImageFunction_h
#define __itkBlur3DImageFunction_h

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
  itkStaticConstMacro(ImageDimension, unsigned int,
    ::itk::GetImageDimension< InputImageType >::ImageDimension ); 

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
  virtual double EvaluateAtContinuousIndex( const ContinuousIndexType &
    index ) const ;

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

  typedef std::list< double >                        KernelWeightsListType;
  typedef std::list< typename InputImageType::IndexType >  
                                                     KernelXListType;

  bool                    m_Debug;

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

  IndexType               m_ImageIndexMin;
  IndexType               m_ImageIndexMax;
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlur3DImageFunction.txx"
#endif

#endif
