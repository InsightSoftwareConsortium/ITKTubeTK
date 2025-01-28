/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#ifndef __itktubePointBasedSpatialObjectTransformFilter_hxx
#define __itktubePointBasedSpatialObjectTransformFilter_hxx


#include <itkSpatialObject.h>
#include <itkSpatialObjectFactory.h>

namespace itk
{

namespace tube
{

template <class TTransformType, unsigned int TDimension>
PointBasedSpatialObjectTransformFilter<TTransformType, TDimension>::PointBasedSpatialObjectTransformFilter(void)
{
  m_OutputObjectToParentTransform = 0;
  m_Transform = 0;

  SpatialObjectFactoryBase::RegisterDefaultSpatialObjects();
  SpatialObjectFactory<SpatialObject<TDimension>>::RegisterSpatialObject();
  SpatialObjectFactory<TubeType>::RegisterSpatialObject();
  SpatialObjectFactory<SurfaceType>::RegisterSpatialObject();
  SpatialObjectFactory<LineType>::RegisterSpatialObject();
  SpatialObjectFactory<DTITubeType>::RegisterSpatialObject();
  SpatialObjectFactory<ContourType>::RegisterSpatialObject();
}

/**
 * Apply the transformation to the tube
 */
template <class TTransformType, unsigned int TDimension>
void
PointBasedSpatialObjectTransformFilter<TTransformType, TDimension>::GenerateData(void)
{
  typename SpatialObject<TDimension>::ConstPointer inputSO = this->GetInput();
  typename SpatialObject<TDimension>::Pointer      outputSO = this->GetOutput();

  Transform(inputSO, outputSO);

  typedef typename SpatialObject<TDimension>::ChildrenConstListType ChildrenConstListType;
  ChildrenConstListType *                                           children = this->GetInput()->GetConstChildren();
  typename ChildrenConstListType::iterator                          it = children->begin();
  while (it != children->end())
  {
    typename SpatialObject<TDimension>::Pointer tmpSO = (*it)->Clone();
    this->UpdateLevel(*it, outputSO);
    ++it;
  }
  delete children;

  this->GraftOutput(outputSO);
}

/**
 * Apply the transformation to the tube
 */
template <class TTransformType, unsigned int TDimension>
void
PointBasedSpatialObjectTransformFilter<TTransformType, TDimension>::UpdateLevel(
  const SpatialObject<TDimension> * inputSO,
  SpatialObject<TDimension> *       parentSO)
{
  typename SpatialObject<TDimension>::Pointer outputSO = inputSO->Clone();

  Transform(inputSO, outputSO);

  parentSO->AddChild(outputSO);

  typedef typename SpatialObject<TDimension>::ChildrenListType ChildrenListType;
  ChildrenListType *                                           children = inputSO->GetChildren();
  typename ChildrenListType::const_iterator                    it = children->begin();
  while (it != children->end())
  {
    this->UpdateLevel(*it, outputSO);
    ++it;
  }
  delete children;
}

template <class TTransformType, unsigned int TDimension>
bool
PointBasedSpatialObjectTransformFilter<TTransformType, TDimension>::Transform(const SpatialObject<TDimension> * inputSO,
                                                                              SpatialObject<TDimension> * outputSO)
{

  const TubeType * inputSOAsTube = dynamic_cast<const TubeSpatialObject<TDimension> *>(inputSO);
  if (inputSOAsTube != nullptr)
  {
    TubeSpatialObject<TDimension> * outputSOAsTube = dynamic_cast<TubeSpatialObject<TDimension> *>(outputSO);
    this->TransformTube(inputSOAsTube, outputSOAsTube);
    return true;
  }

  const SurfaceType * inputSOAsSurface = dynamic_cast<const SurfaceSpatialObject<TDimension> *>(inputSO);
  if (inputSOAsSurface != nullptr)
  {
    SurfaceSpatialObject<TDimension> * outputSOAsSurface = dynamic_cast<SurfaceSpatialObject<TDimension> *>(outputSO);
    this->TransformSurface(inputSOAsSurface, outputSOAsSurface);
    return true;
  }

  const PointBasedType * inputSOAsPointBased = dynamic_cast<const PointBasedSpatialObject<TDimension> *>(inputSO);
  if (inputSOAsPointBased != nullptr)
  {
    PointBasedSpatialObject<TDimension> * outputSOAsPointBased =
      dynamic_cast<PointBasedSpatialObject<TDimension> *>(outputSO);
    this->TransformPointBased(inputSOAsPointBased, outputSOAsPointBased);
    return true;
  }

  return false;
}

template <class TTransformType, unsigned int TDimension>
void
PointBasedSpatialObjectTransformFilter<TTransformType, TDimension>::TransformPointBased(
  const PointBasedSpatialObject<TDimension> * inputSO,
  PointBasedSpatialObject<TDimension> *       outputSO)
{
  Point<double, TDimension> objectPoint;
  Point<double, TDimension> transformedObjectPoint;

  outputSO->CopyInformation(inputSO);
  outputSO->Clear();

  if (m_OutputObjectToParentTransform.IsNotNull())
  {
    typename SpatialObjectTransformType::Pointer tfm = SpatialObjectTransformType::New();
    tfm->SetIdentity();
    tfm->SetMatrix(m_OutputObjectToParentTransform->GetMatrix());
    tfm->SetOffset(m_OutputObjectToParentTransform->GetOffset());
    outputSO->SetObjectToParentTransform(tfm);
  }

  typedef typename PointBasedType::SpatialObjectPointListType SpatialObjectPointListType;
  SpatialObjectPointListType                                  pointBasedPointList = inputSO->GetPoints();
  typename SpatialObjectPointListType::const_iterator         pointBasedPointIterator = pointBasedPointList.begin();

  while (pointBasedPointIterator != pointBasedPointList.end())
  {
    objectPoint = (*pointBasedPointIterator).GetPositionInObjectSpace();

    transformedObjectPoint = m_Transform->TransformPoint(objectPoint);

    SpatialObjectPoint<TDimension> pnt;

    pnt.SetId(pointBasedPointIterator->GetId());
    pnt.SetColor(pointBasedPointIterator->GetColor());

    pnt.SetPositionInObjectSpace(transformedObjectPoint);

    outputSO->AddPoint(pnt);

    ++pointBasedPointIterator;
  }
}

template <class TTransformType, unsigned int TDimension>
void
PointBasedSpatialObjectTransformFilter<TTransformType, TDimension>::TransformTube(
  const TubeSpatialObject<TDimension> * inputSO,
  TubeSpatialObject<TDimension> *       outputSO)
{
  Point<double, TDimension> objectPoint;
  Point<double, TDimension> transformedObjectPoint;

  outputSO->Clear();
  outputSO->CopyInformation(inputSO);

  if (m_OutputObjectToParentTransform.IsNotNull())
  {
    typename SpatialObjectTransformType::Pointer tfm = SpatialObjectTransformType::New();
    tfm->SetIdentity();
    tfm->SetMatrix(m_OutputObjectToParentTransform->GetMatrix());
    tfm->SetOffset(m_OutputObjectToParentTransform->GetOffset());
    outputSO->SetObjectToParentTransform(tfm);
  }

  typedef typename TubeType::TubePointListType TubePointListType;
  TubePointListType                            tubePointList = inputSO->GetPoints();
  typename TubePointListType::const_iterator   tubePointIterator = tubePointList.begin();

  while (tubePointIterator != tubePointList.end())
  {
    objectPoint = (*tubePointIterator).GetPositionInObjectSpace();

    transformedObjectPoint = m_Transform->TransformPoint(objectPoint);

    TubeSpatialObjectPoint<TDimension> pnt;

    pnt.SetId(tubePointIterator->GetId());
    pnt.SetColor(tubePointIterator->GetColor());

    pnt.SetSpatialObject(outputSO);

    pnt.SetPositionInObjectSpace(transformedObjectPoint);

    // get both normals
    typename TubeType::CovariantVectorType n1 = tubePointIterator->GetNormal1InObjectSpace();
    typename TubeType::CovariantVectorType n2 = tubePointIterator->GetNormal2InObjectSpace();

    // only try transformation of normals if both are non-zero
    if (!n1.GetVnlVector().is_zero() && !n2.GetVnlVector().is_zero())
    {
      n1 = m_Transform->TransformCovariantVector(n1, objectPoint);
      n2 = m_Transform->TransformCovariantVector(n2, objectPoint);
      n1.Normalize();
      n2.Normalize();
      pnt.SetNormal1InObjectSpace(n1);
      pnt.SetNormal2InObjectSpace(n2);
    }

    typename TubeType::VectorType tang = tubePointIterator->GetTangentInObjectSpace();
    if (!tang.GetVnlVector().is_zero())
    {
      tang = m_Transform->TransformVector(tang, objectPoint);
      tang.Normalize();
      pnt.SetTangentInObjectSpace(tang);
    }

    typename TubeType::VectorType radi;
    for (unsigned int i = 0; i < TDimension; ++i)
    {
      radi[i] = tubePointIterator->GetRadiusInObjectSpace();
    }
    radi = m_Transform->TransformVector(radi, objectPoint);
    pnt.SetRadiusInObjectSpace(radi[0]);

    pnt.SetMedialness((*tubePointIterator).GetMedialness());
    pnt.SetRidgeness((*tubePointIterator).GetRidgeness());
    pnt.SetBranchness((*tubePointIterator).GetBranchness());

    outputSO->AddPoint(pnt);

    ++tubePointIterator;
  }
}


template <class TTransformType, unsigned int TDimension>
void
PointBasedSpatialObjectTransformFilter<TTransformType, TDimension>::TransformSurface(
  const SurfaceSpatialObject<TDimension> * inputSO,
  SurfaceSpatialObject<TDimension> *       outputSO)
{
  // We make the copy and sub-sample if it is a tube.
  Point<double, TDimension> objectPoint;
  Point<double, TDimension> transformedObjectPoint;

  outputSO->CopyInformation(inputSO);
  outputSO->Clear();

  if (m_OutputObjectToParentTransform.IsNotNull())
  {
    typename SpatialObjectTransformType::Pointer tfm = SpatialObjectTransformType::New();
    tfm->SetIdentity();
    tfm->SetMatrix(m_OutputObjectToParentTransform->GetMatrix());
    tfm->SetOffset(m_OutputObjectToParentTransform->GetOffset());
    outputSO->SetObjectToParentTransform(tfm);
  }

  typedef typename SurfaceType::SurfacePointListType SurfacePointListType;
  SurfacePointListType                               surfacePointList = inputSO->GetPoints();
  typename SurfacePointListType::const_iterator      surfacePointIterator = surfacePointList.begin();

  while (surfacePointIterator != surfacePointList.end())
  {
    objectPoint = (*surfacePointIterator).GetPositionInObjectSpace();

    transformedObjectPoint = m_Transform->TransformPoint(objectPoint);

    SurfaceSpatialObjectPoint<TDimension> pnt;

    pnt.SetId(surfacePointIterator->GetId());
    pnt.SetColor(surfacePointIterator->GetColor());

    pnt.SetPositionInObjectSpace(transformedObjectPoint);

    // get normals
    typename SurfaceType::CovariantVectorType n1 = surfacePointIterator->GetNormalInObjectSpace();

    // only try transformation of normals if both are non-zero
    if (!n1.GetVnlVector().is_zero())
    {
      n1 = m_Transform->TransformCovariantVector(n1, objectPoint);
      n1.Normalize();
      pnt.SetNormalInObjectSpace(n1);
    }

    outputSO->AddPoint(pnt);

    ++surfacePointIterator;
  }
}

template <class TTransformType, unsigned int TDimension>
void
PointBasedSpatialObjectTransformFilter<TTransformType, TDimension>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Transformation: " << m_Transform << std::endl;
  os << indent << "OutputObjectToParent Transform: " << m_OutputObjectToParentTransform << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubePointBasedSpatialObjectTransformFilter_hxx )
