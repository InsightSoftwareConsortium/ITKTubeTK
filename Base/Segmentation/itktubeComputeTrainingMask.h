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

#ifndef __itktubeComputeTrainingMask_h
#define __itktubeComputeTrainingMask_h

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
#include <itkDivideImageFilter.h>

namespace itk
{

namespace tube
{

/**
 * This class returns expert vessel and not vessel mask.
 *
 * \sa ComputeTrainingMask
 */

template< class TInputImage >
class ComputeTrainingMask:
        public ImageToImageFilter< TInputImage,
                                   itk::Image<short,TInputImage::ImageDimension> >
{
public:
  typedef ComputeTrainingMask                             Self;
  typedef ImageToImageFilter<TInputImage,TInputImage>     Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;
  typedef TInputImage                                     ImageType;
  typedef itk::Image< short, ImageType::ImageDimension >  ImageTypeShort;

  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  itkNewMacro( Self );
  const ImageTypeShort* GetNotVesselMask();
  itkSetMacro(Gap,double);
  itkSetMacro(NotVesselWidth,double);

protected:
  ComputeTrainingMask();
  virtual ~ComputeTrainingMask();
  virtual void GenerateData();
  void PrintSelf( std::ostream & os, Indent indent ) const;

private:
  typedef itk::BinaryBallStructuringElement< short, ImageType::ImageDimension > BallType;

  typedef itk::DilateObjectMorphologyImageFilter< ImageType, ImageType, BallType >
  DilateFilterType;
  typedef itk::BinaryThinningImageFilter< ImageType, ImageType >
  BinaryThinningFilterType;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType >
  ThresholdFilterType;
  typedef itk::SubtractImageFilter<ImageType,ImageType,ImageType>
  SubstractFilterType;
  typedef itk::MultiplyImageFilter<ImageType,ImageType,ImageType>
  MultiplyFilterType;
  typedef itk::DivideImageFilter<ImageType,ImageType,ImageType>
  DivideFilterType;
  typedef itk::AddImageFilter<ImageType,ImageType,ImageType>
  AddFilterType;
  typedef itk::CastImageFilter< ImageType, ImageTypeShort >
  CastFilterType;

  ComputeTrainingMask( const Self& );
  void operator=( const Self& );
  void ApplyDilateMorphologyFilter( typename ImageType::Pointer &input );

  typename AddFilterType::Pointer             m_Add;
  typename ThresholdFilterType::Pointer       m_Threshold;
  typename BinaryThinningFilterType::Pointer  m_BinaryThinning;
  typename DilateFilterType::Pointer          m_Dilate;
  typename SubstractFilterType::Pointer       m_Substract;
  typename MultiplyFilterType::Pointer        m_MultiplyCenterLine;
  typename DivideFilterType::Pointer          m_DivideImage;
  typename CastFilterType::Pointer            m_Cast;
  typename CastFilterType::Pointer            m_CastNotVessel;

  BallType m_Ball;
  double   m_Gap;
  double   m_NotVesselWidth;
};

}//end of tube namespace
}//end of itk namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeComputeTrainingMask.hxx"
#endif

#endif
