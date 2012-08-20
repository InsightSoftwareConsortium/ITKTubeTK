/*=========================================================================

Library:   TubeTK

Copyright 2012 Kitware Inc. 28 Corporate Drive,
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
#ifndef __itkAngleOfIncidenceImageFilter_h
#define __itkAngleOfIncidenceImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkSymmetricEigenAnalysisImageFilter.h"
#include "itkSymmetricEigenVectorAnalysisImageFilter.h"
#include "itkVectorImage.h"
#include "itkHessianRecursiveGaussianImageFilter.h"

namespace itk
{
/** \class AngleOfIncidenceImageFilter
 * \brief Computes angle of incidence
 *
 * The angle of incidence is defined as the angle between the beam
 * direction at a organ boundary and the normal to the boundary.
 *
 * \ingroup ImageToImageFilter
 */
template< class TInputImage, class TOutputImage >
class ITK_EXPORT AngleOfIncidenceImageFilter:public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef AngleOfIncidenceImageFilter                         Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage >     Superclass;
  typedef SmartPointer< Self >                                Pointer;
  typedef SmartPointer< const Self >                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(AngleOfIncidenceImageFilter, ImageToImageFilter);

  itkStaticConstMacro(Dimension, unsigned int, 3);

  /** Some convenient typedefs for input image */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::ConstPointer InputImagePointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;
  typedef typename InputImageType::PixelType    InputImagePixelType;


  /** typedef for the origin type */
  typedef Point< double, 3>     VectorType;

  /* typedef for the output image */
  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::Pointer    OutputImagePointer;
  typedef typename OutputImageType::RegionType OutputImageRegionType;
  typedef typename OutputImageType::PixelType  OutputImagePixelType;

  /** typedef to generate surface normal vector using eigen analysis */
  typedef typename itk::HessianRecursiveGaussianImageFilter< InputImageType >  HessianFilterType;

  typedef itk::SymmetricSecondRankTensor< double, 3 >         SymmetricSecondRankTensorType;
  typedef itk::Image< SymmetricSecondRankTensorType, 3>       SymmetricSecondRankTensorImageType;
  typedef  itk::Matrix< double, 3, 3>                         EigenVectorMatrixType;
  typedef  itk::Image< EigenVectorMatrixType, 3>              EigenVectorMatrixImageType;
  typedef  itk::FixedArray< double, 3>                        EigenValueArrayType;
  typedef  itk::Image< EigenValueArrayType, 3>                EigenValueImageType;

  typedef itk::VectorImage< double, 3 >    EigenVectorImageType;

  typedef itk::SymmetricEigenAnalysisImageFilter
        <SymmetricSecondRankTensorImageType, EigenValueImageType> EigenValueAnalysisFilterType;

  typedef itk::SymmetricEigenVectorAnalysisImageFilter
        <SymmetricSecondRankTensorImageType, EigenValueImageType,
         EigenVectorMatrixImageType>                                    EigenVectorAnalysisFilterType;

  /** Set/Get Ultrasound origin vector */
  itkSetMacro(UltrasoundOrigin, VectorType);
  itkGetConstMacro(UltrasoundOrigin, VectorType);

  /* Generate Data */
  void GenerateData(void);

protected:
  AngleOfIncidenceImageFilter();
  virtual ~AngleOfIncidenceImageFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  void ComputeNormalVectorImage();
private:
  AngleOfIncidenceImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);          //purposely not implemented

  /* Ultasound origin*/
  VectorType m_UltrasoundOrigin;

  /* Hessian analysis filter */
  typename HessianFilterType::Pointer m_HessianFilter;

  /* Eigen value analysis filter */
  EigenValueAnalysisFilterType::Pointer m_EigenValueAnalysisFilter;

  /* Eigen vector analysis filter */
  EigenVectorAnalysisFilterType::Pointer m_EigenVectorAnalysisFilter;

  // Primary eigen vector image
  EigenVectorImageType::Pointer m_PrimaryEigenVectorImage;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAngleOfIncidenceImageFilter.txx"
#endif

#endif
