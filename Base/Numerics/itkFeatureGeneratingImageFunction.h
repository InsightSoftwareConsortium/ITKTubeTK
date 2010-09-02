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
#ifndef __itkFeatureGeneratingImageFunction_h
#define __itkFeatureGeneratingImageFunction_h

#include <vector>
#include <map>

#include "itkImageFunction.h"
#include "itkOrientedImage.h"

namespace itk
{

/** \class FeatureGeneratingImageFunction
 *
 */
template<class TInputImage, class TCoordRep = float>
class ITK_EXPORT FeatureGeneratingImageFunction :
  public ImageFunction< TInputImage, std::vector<double>, TCoordRep >
{
public:

  /** Class typedefs **/
  typedef FeatureGeneratingImageFunction               Self;
  typedef ImageFunction<TInputImage,std::vector<double>,TCoordRep>  
                                                       Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;
  typedef typename Superclass::InputImageType          InputImageType;
  typedef typename TInputImage::PixelType              PixelType;
  typedef typename Superclass::PointType               PointType;
  typedef typename Superclass::IndexType               IndexType;
  typedef typename Superclass::ContinuousIndexType     ContinuousIndexType;
  typedef std::vector<double>                          OutputType;
  typedef std::vector<std::string>                     FeatureListType;

  /** Run-time type information (and related methods). */
  itkTypeMacro( FeatureGeneratingImageFunction, ImageFunction );

  /** Standard New Macro. */
  itkNewMacro( Self );

  /** Constant for fetching the dimensions of the image. **/
  itkStaticConstMacro( ImageDimension, unsigned int,
                       Superclass::ImageDimension );

  /** Override the Set for the InputImage */
  virtual void SetInputImage( const InputImageType * ptr );

  /** Get the output vector at a given point. */
  virtual OutputType Evaluate( const PointType& point ) const
    {
    IndexType index;
    this->ConvertPointToNearestIndex( point, index );
    return ( this->EvaluateAtIndex( index ) );
    }

  /** Get the output vector at a given continuous index. */
  virtual OutputType EvaluateAtContinuousIndex( 
    const ContinuousIndexType & index ) const
    {
    IndexType nindex;

    this->ConvertContinuousIndexToNearestIndex( index, nindex );
    return this->EvaluateAtIndex( nindex );
    }

  /** Get the feature vector at an index for a given point **/
  virtual OutputType EvaluateAtIndex( const IndexType & index ) const;

  /** Get the output vector at a given point. */
  virtual std::string EvaluateToString( const PointType& point ) const
    {
    IndexType index;
    this->ConvertPointToNearestIndex( point, index );
    return ( this->EvaluateToStringAtIndex( index ) );
    }

  /** Get the output vector at a given continuous index. */
  virtual std::string EvaluateToStringAtContinuousIndex( 
    const ContinuousIndexType & index ) const
    {
    IndexType nindex;

    this->ConvertContinuousIndexToNearestIndex( index, nindex );
    return this->EvaluateToStringAtIndex( nindex );
    }

  /** Get the feature vector at an index for a given point **/
  std::string EvaluateToStringAtIndex( const IndexType & index ) const;

  virtual const FeatureListType& GetFeatureLabels() const;

protected:
  
  /** Default constructor */
  FeatureGeneratingImageFunction();

  /** Default destructor */
  ~FeatureGeneratingImageFunction() {}

  /** Printself function for introspection. **/
  void PrintSelf( std::ostream& os, Indent indent ) const;

  FeatureListType m_Features;

private:
  FeatureGeneratingImageFunction( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

}; // End class FeatureGeneratingImageFunction

}// end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkFeatureGeneratingImageFunction.txx"
#endif

#endif
