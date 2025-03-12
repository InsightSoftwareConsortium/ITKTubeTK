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

#ifndef itktubeBasisFeatureVectorGenerator_h
#define itktubeBasisFeatureVectorGenerator_h

#include "itktubeFeatureVectorGenerator.h"

#include <itkImage.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include <vector>

namespace itk
{

namespace tube
{

template <class TImage, class TLabelMap>
class BasisFeatureVectorGenerator : public FeatureVectorGenerator<TImage>
{
public:
  using Self = BasisFeatureVectorGenerator;
  using Superclass = FeatureVectorGenerator<TImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  itkTypeMacro(BasisFeatureVectorGenerator, FeatureVectorGenerator);

  itkNewMacro(Self);

  //
  using ImageType = typename Superclass::ImageType;

  using LabelMapType = TLabelMap;

  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);

  using IndexType = typename Superclass::IndexType;

  using FeatureValueType = typename Superclass::FeatureValueType;
  using FeatureVectorType = typename Superclass::FeatureVectorType;
  using FeatureImageType = typename Superclass::FeatureImageType;

  using FeatureVectorGeneratorType = FeatureVectorGenerator<TImage>;

  using ObjectIdType = typename TLabelMap::PixelType;
  using ObjectIdListType = std::vector<ObjectIdType>;

  using ValueType = typename Superclass::ValueType;
  using ValueListType = typename Superclass::ValueListType;

  using VectorType = vnl_vector<ValueType>;
  using VectorListType = std::vector<VectorType>;

  using MatrixType = vnl_matrix<ValueType>;
  using MatrixListType = std::vector<MatrixType>;

  void
  SetInputFeatureVectorGenerator(FeatureVectorGeneratorType * fGen);
  typename FeatureVectorGenerator<TImage>::Pointer
  GetInputFeatureVectorGenerator(void);

  itkSetObjectMacro(LabelMap, LabelMapType);
  itkGetModifiableObjectMacro(LabelMap, LabelMapType);

  void
  SetObjectId(ObjectIdType objectId);
  void
  AddObjectId(ObjectIdType objectId);
  ObjectIdType
  GetObjectId(unsigned int num = 0) const;
  unsigned int
  GetNumberOfObjectIds(void) const;

  ValueType
  GetObjectMean(ObjectIdType objectId) const;
  void
  SetObjectMean(ObjectIdType objectId, ValueType val);
  MatrixType
  GetObjectCovariance(ObjectIdType objectId) const;
  void
  SetObjectCovariance(ObjectIdType objectId, MatrixType val);

  ValueType
  GetGlobalMean(void) const;
  void
  SetGlobalMean(ValueType val);
  MatrixType
  GetGlobalCovariance(void) const;
  void
  SetGlobalCovariance(MatrixType val);

  void
  SetNumberOfPCABasisToUseAsFeatures(unsigned int numBasisUsed);
  unsigned int
  GetNumberOfPCABasisToUseAsFeatures(void) const;
  void
  SetNumberOfLDABasisToUseAsFeatures(unsigned int numBasisUsed);
  unsigned int
  GetNumberOfLDABasisToUseAsFeatures(void) const;

  double
  GetBasisValue(unsigned int basisNum) const;
  VectorType
  GetBasisVector(unsigned int basisNum) const;
  MatrixType
  GetBasisMatrix(void) const;
  VectorType
  GetBasisValues(void) const;

  void
  SetBasisValue(unsigned int basisNum, double value);
  void
  SetBasisVector(unsigned int basisNum, const VectorType & vec);
  void
  SetBasisMatrix(const MatrixType & mat);
  void
  SetBasisValues(const VectorType & values);

  virtual typename FeatureImageType::Pointer
  GetFeatureImage(unsigned int fNum) const override;

  void
  SetInputWhitenMeans(const ValueListType & means);
  const ValueListType &
  GetInputWhitenMeans(void) const;
  void
  SetInputWhitenStdDevs(const ValueListType & stdDevs);
  const ValueListType &
  GetInputWhitenStdDevs(void) const;

  void
  SetOutputWhitenMeans(const ValueListType & means);
  const ValueListType &
  GetOutputWhitenMeans(void) const;
  void
  SetOutputWhitenStdDevs(const ValueListType & stdDevs);
  const ValueListType &
  GetOutputWhitenStdDevs(void) const;

  virtual void
  Update(void) override;

  virtual unsigned int
  GetNumberOfFeatures(void) const override;

  virtual FeatureVectorType
  GetFeatureVector(const IndexType & indx) const override;

  virtual FeatureValueType
  GetFeatureVectorValue(const IndexType & indx, unsigned int fNum) const override;

protected:
  BasisFeatureVectorGenerator(void);
  virtual ~BasisFeatureVectorGenerator(void);

  virtual void
  UpdateWhitenStatistics(void) override;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  // Purposely not implemented
  BasisFeatureVectorGenerator(const Self &);
  void
  operator=(const Self &); // Purposely not implemented

  //  Data
  typename FeatureVectorGeneratorType::Pointer m_InputFeatureVectorGenerator;

  typename LabelMapType::Pointer m_LabelMap;

  ObjectIdListType m_ObjectIdList;
  VectorListType   m_ObjectMeanList;
  MatrixListType   m_ObjectCovarianceList;

  VectorType m_GlobalMean;
  MatrixType m_GlobalCovariance;

  unsigned int m_NumberOfPCABasisToUseAsFeatures;
  unsigned int m_NumberOfLDABasisToUseAsFeatures;

  MatrixType m_BasisMatrix;
  VectorType m_BasisValues;
}; // End class BasisFeatureVectorGenerator

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeBasisFeatureVectorGenerator.hxx"
#endif

#endif // End !defined( __itktubeBasisFeatureVectorGenerator_h )
