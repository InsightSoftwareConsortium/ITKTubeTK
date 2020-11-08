/*=========================================================================
   Library:   TubeTK

   Copyright Kitware Inc.

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

#ifndef __itktubeComputeTrainingMaskFilter_h
#define __itktubeComputeTrainingMaskFilter_h

#include <itkImageToImageFilter.h>

#include <itkBinaryThinningImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkDilateObjectMorphologyImageFilter.h>
#include <itkErodeObjectMorphologyImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>

namespace itk
{

namespace tube
{

/**
 * This class returns expert vessel and not vessel mask.
 *
 * \sa ComputeTrainingMaskFilter
 */

template< class TInputImage, class TLabelMap=Image<
  typename TInputImage::PixelType, TInputImage::ImageDimension> >
class ComputeTrainingMaskFilter:
  public ImageToImageFilter< TInputImage,
    TLabelMap >
{
public:
  typedef ComputeTrainingMaskFilter                       Self;
  typedef ImageToImageFilter<TInputImage, TLabelMap>      Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;
  typedef TInputImage                                     ImageType;
  typedef TLabelMap                                       LabelMapType;

  itkStaticConstMacro( InputImageDimension, unsigned int,
                      TInputImage::ImageDimension );

  itkNewMacro( Self );
  const LabelMapType * GetObjectMask();
  const LabelMapType * GetNotObjectMask();
  itkSetMacro( Gap, double );
  itkSetMacro( ObjectWidth, double );
  itkSetMacro( NotObjectWidth, double );
  itkGetMacro( Gap, double );
  itkGetMacro( ObjectWidth, double );
  itkGetMacro( NotObjectWidth, double );

protected:
  ComputeTrainingMaskFilter();
  virtual ~ComputeTrainingMaskFilter();

  virtual void GenerateData() override;

  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:
  typedef itk::BinaryBallStructuringElement< short,
    ImageType::ImageDimension > BallType;
  typedef itk::DilateObjectMorphologyImageFilter< ImageType, ImageType,
    BallType >                  DilateFilterType;
  typedef itk::BinaryThinningImageFilter< ImageType, ImageType >
                                BinaryThinningFilterType;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType >
                                ThresholdFilterType;
  typedef itk::SubtractImageFilter<ImageType, ImageType, ImageType>
                                SubtractFilterType;
  typedef itk::MultiplyImageFilter<ImageType, ImageType, ImageType>
                                MultiplyFilterType;
  typedef itk::AddImageFilter<ImageType, ImageType, ImageType>
                                AddFilterType;
  typedef itk::CastImageFilter< ImageType, LabelMapType >
                                CastFilterType;

  ComputeTrainingMaskFilter( const Self& );
  void operator=( const Self& );
  void ApplyDilateMorphologyFilter( typename ImageType::Pointer &input,
    int size=1 );

  typename AddFilterType::Pointer             m_Add;
  typename MultiplyFilterType::Pointer        m_Multiply;
  typename ThresholdFilterType::Pointer       m_Threshold;
  typename BinaryThinningFilterType::Pointer  m_BinaryThinning;
  typename DilateFilterType::Pointer          m_Dilate;
  typename SubtractFilterType::Pointer        m_Subtract;
  typename MultiplyFilterType::Pointer        m_MultiplyCenterLine;
  typename MultiplyFilterType::Pointer        m_MultiplyOutside;
  typename CastFilterType::Pointer            m_Cast;
  typename CastFilterType::Pointer            m_CastObject;
  typename CastFilterType::Pointer            m_CastNotObject;

  BallType m_Ball;
  double   m_Gap;
  double   m_ObjectWidth;
  double   m_NotObjectWidth;
};

}//end of tube namespace
}//end of itk namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeComputeTrainingMaskFilter.hxx"
#endif

#endif
