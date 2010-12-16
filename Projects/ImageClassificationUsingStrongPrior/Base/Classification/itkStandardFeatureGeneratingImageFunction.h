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
#ifndef __itkStandardFeatureGeneratingImageFunction_h
#define __itkStandardFeatureGeneratingImageFunction_h

#include <vector>

#include "itkImageFunction.h"
#include "itkOrientedImage.h"

#include "itkTubeRidgeExtractor.h"
#include "itkJointHistogramImageFunction.h"

#include "itkFeatureGeneratingImageFunction.h"

namespace itk
{

/** \class StandardFeatureGeneratingImageFunction
 *
 */
template<class TInputImage, class TCoordRep = float>
class ITK_EXPORT StandardFeatureGeneratingImageFunction :
  public FeatureGeneratingImageFunction< TInputImage, TCoordRep >
{
public:

  /** Class typedefs **/
  typedef StandardFeatureGeneratingImageFunction       Self;
  typedef FeatureGeneratingImageFunction<TInputImage,TCoordRep>  
                                                       Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;
  typedef typename Superclass::InputImageType          InputImageType;
  typedef typename TInputImage::PixelType              PixelType;
  typedef typename Superclass::PointType               PointType;
  typedef typename Superclass::IndexType               IndexType;
  typedef typename Superclass::ContinuousIndexType     ContinuousIndexType;
  typedef typename Superclass::OutputType              OutputType;
  typedef typename Superclass::FeatureListType         FeatureListType;
  typedef itk::tube::RidgeExtractor< InputImageType >  CalculatorType;
  typedef itk::JointHistogramImageFunction<InputImageType>
                                                       HistCalcType;
  typedef typename HistCalcType::HistogramType         HistogramType;

  /** Run-time type information (and related methods). */
  itkTypeMacro( StandardFeatureGeneratingImageFunction, 
    FeatureGeneratingImageFunction );

  /** Standard New Macro. */
  itkNewMacro( Self );

  /** Constant for fetching the dimensions of the image. **/
  itkStaticConstMacro( ImageDimension, unsigned int,
    Superclass::ImageDimension );
  
  void SetPriorImage( typename InputImageType::Pointer prior )
    {
    m_Prior = prior;
    };

  void SetMeanAddHistogram( typename HistogramType::Pointer addMeanHist )
    {
    m_AddMeanHist = addMeanHist;
    };

  void SetStdevAddHistogram( typename HistogramType::Pointer addStdevHist )
    {
    m_AddStdevHist = addStdevHist;
    };

  void SetMeanSubHistogram( typename HistogramType::Pointer subMeanHist )
    {
    m_SubMeanHist = subMeanHist;
    };

  void SetStdevSubHistogram( typename HistogramType::Pointer subStdevHist )
    {
    m_SubStdevHist = subStdevHist;
    };

  void SetMeanNormHistogram( typename HistogramType::Pointer normMeanHist )
    {
    m_NormMeanHist = normMeanHist;
    };

  void SetStdevNormHistogram( typename HistogramType::Pointer normStdevHist )
    {
    m_NormStdevHist = normStdevHist;
    };

  void SetScale( double scale )
    {
    m_Scale = scale;
    m_SigmaMedium = (m_Scale/2)*0.6667;
    m_SigmaSmall = 0.6667*m_SigmaMedium;
    m_SigmaLarge = 1.3333*m_SigmaMedium;
    };

  /** Prepare the filter to go **/
  void PrepFilter();
  
  /** Get the feature vector at an index for a given point **/
  virtual OutputType EvaluateAtIndex( const IndexType & index ) const;

protected:
  
  /** Default constructor */
  StandardFeatureGeneratingImageFunction();

  /** Default destructor */
  ~StandardFeatureGeneratingImageFunction() {};

  typename InputImageType::Pointer m_Prior;
  typename HistogramType::Pointer  m_AddMeanHist;
  typename HistogramType::Pointer  m_AddStdevHist;
  typename HistogramType::Pointer  m_SubMeanHist;
  typename HistogramType::Pointer  m_SubStdevHist;
  typename HistogramType::Pointer  m_NormMeanHist;
  typename HistogramType::Pointer  m_NormStdevHist;

  typename CalculatorType::Pointer m_InputCalcSmall;
  typename CalculatorType::Pointer m_InputCalcMedium;
  typename CalculatorType::Pointer m_InputCalcLarge;
  typename CalculatorType::Pointer m_PriorCalcSmall;
  typename CalculatorType::Pointer m_PriorCalcMedium;
  typename CalculatorType::Pointer m_PriorCalcLarge;

  typename HistCalcType::Pointer m_AddJHCalc;
  typename HistCalcType::Pointer m_SubJHCalc;
  typename HistCalcType::Pointer m_NomJHCalc;

  double m_InputDataMin;
  double m_InputDataMax;
  double m_PriorDataMin;
  double m_PriorDataMax;
  double m_Scale;

  PixelType m_SigmaMedium;
  PixelType m_SigmaSmall;
  PixelType m_SigmaLarge;

private:
  StandardFeatureGeneratingImageFunction( const Self& ); //pni
  void operator=( const Self& ); //purposely not implemented

}; // End class StandardFeatureGeneratingImageFunction

}// end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkStandardFeatureGeneratingImageFunction.txx"
#endif

#endif
