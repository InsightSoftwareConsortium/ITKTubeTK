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

#ifndef itktubeCVTImageFilter_h
#define itktubeCVTImageFilter_h

#include <itkContinuousIndex.h>
#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkIndex.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkProcessObject.h>
#include <itkVariableSizeMatrix.h>

#include <vector>

namespace itk
{

namespace tube
{

template <class TInputImage, class TOutputImage = TInputImage>
class CVTImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class type alias. */
  using Self = CVTImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkOverrideGetNameOfClassMacro(CVTImageFilter);

  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  using InputImageType = TInputImage;
  using InputPixelType = typename InputImageType::PixelType;

  using IndexType = typename InputImageType::IndexType;
  using ContinuousIndexType = ContinuousIndex<double, itkGetStaticConstMacro(ImageDimension)>;

  using OutputImageType = TOutputImage;
  using OutputPixelType = typename OutputImageType::PixelType;

  using RegionType = typename InputImageType::RegionType;
  using SizeType = typename RegionType::SizeType;

  using PointArrayType = std::vector<ContinuousIndexType>;

  using SamplingMethodEnum = enum { CVT_GRID, CVT_RANDOM, CVT_USER };

  /** */
  itkGetMacro(NumberOfCentroids, unsigned int);
  itkSetMacro(NumberOfCentroids, unsigned int);

  /** */
  itkGetMacro(InitialSamplingMethod, SamplingMethodEnum);
  itkSetMacro(InitialSamplingMethod, SamplingMethodEnum);

  /** */
  itkGetMacro(NumberOfSamples, unsigned int);
  itkSetMacro(NumberOfSamples, unsigned int);

  /** */
  itkGetMacro(NumberOfIterations, unsigned int);
  itkSetMacro(NumberOfIterations, unsigned int);

  /** */
  itkGetMacro(NumberOfIterationsPerBatch, unsigned int);
  itkSetMacro(NumberOfIterationsPerBatch, unsigned int);

  /** */
  itkGetMacro(NumberOfSamplesPerBatch, unsigned int);
  itkSetMacro(NumberOfSamplesPerBatch, unsigned int);

  /** */
  itkGetMacro(BatchSamplingMethod, SamplingMethodEnum);
  itkSetMacro(BatchSamplingMethod, SamplingMethodEnum);

  /** */
  itkGetMacro(Centroids, PointArrayType);
  void
  SetCentroids(const PointArrayType & centroids);

  /** */
  itkGetMacro(Seed, long int);
  itkSetMacro(Seed, long int);

  itkGetMacro(AdjacencyMatrix, VariableSizeMatrix<double>);

protected:
  CVTImageFilter(void);
  ~CVTImageFilter(void) {}

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  GenerateInputRequestedRegion(void) override;
  void
  EnlargeOutputRequestedRegion(DataObject * output) override;
  void
  GenerateData(void) override;

  double
  ComputeIteration(double & energyDiff);
  void
  ComputeSample(PointArrayType * sample, unsigned int sampleSize, SamplingMethodEnum samplingMethod);
  void
  ComputeClosest(const PointArrayType & sample, const PointArrayType & centroids, unsigned int * nearest);

  void
  ComputeAdjacencyMatrix(void);

private:
  CVTImageFilter(const Self &);
  void
  operator=(const Self &);

  typename OutputImageType::Pointer m_OutputImage;

  typename InputImageType::ConstPointer m_InputImage;

  unsigned int   m_NumberOfCentroids;
  PointArrayType m_Centroids;

  double   m_InputImageMax;
  SizeType m_InputImageSize;

  long int                                                        m_Seed;
  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer m_RandomGenerator;

  SamplingMethodEnum m_InitialSamplingMethod;

  unsigned int m_NumberOfSamples;
  unsigned int m_NumberOfIterations;

  SamplingMethodEnum m_BatchSamplingMethod;
  unsigned int       m_NumberOfIterationsPerBatch;
  unsigned int       m_NumberOfSamplesPerBatch;

  VariableSizeMatrix<double> m_AdjacencyMatrix;

}; // End class CVTImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeCVTImageFilter.hxx"
#endif

#endif // End !defined( _itkCVTImageFilter_h )
