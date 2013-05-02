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
#ifndef __itkTubeToTubeTransformFilter_txx
#define __itkTubeToTubeTransformFilter_txx

#include "itkTubeToTubeTransformFilter.h"

namespace itk
{

/**
 * Constructor
 */
template <class TTransformType, unsigned int TDimension>
TubeToTubeTransformFilter<TTransformType,TDimension>
::TubeToTubeTransformFilter( void )
{
  m_Output = 0;
  m_CropSize = NULL;
  m_NarrowBandSize=0;
  m_Scale=1;
  m_Crop = false;
  m_Transform = 0;
  m_TransformAsGroup = 0;
  m_Ridgeness = -1;
  m_Medialness = -1;
}


/**
 * Apply the transformation to the tube
 */
template <class TTransformType, unsigned int TDimension>
void
TubeToTubeTransformFilter<TTransformType,TDimension>
::Update( void )
{
  m_Output = GroupType::New();

  // Set the spacing first;
  double* groupspacing = new double[TDimension];
  const double* tubespacing;

  for(unsigned int i=0;i<TDimension;i++)
    {
    groupspacing[i]= this->GetInput()->GetIndexToObjectTransform()
                                ->GetScaleComponent()[i]/m_Scale;
    }

  m_Output->GetIndexToObjectTransform()->SetScaleComponent(groupspacing);

  // Set the transform
  if(!m_Transform)
    {
    m_Transform = TransformType::New();
    }

  Point<double,TDimension> point;
  CovariantVector<double,TDimension> FirstNormal;
  CovariantVector<double,TDimension> SecondNormal;

  typename TubeType::ChildrenListType::iterator TubeIterator;

  typedef typename TubeType::PointListType TubePointListType;


  typename TubeType::ChildrenListPointer TubeList =
    this->GetInput()->GetChildren(99999);
  for(TubeIterator = TubeList->begin();
      TubeIterator != TubeList->end();
      TubeIterator++)
    {

    if(!strcmp((*TubeIterator)->GetTypeName(),"VesselTubeSpatialObject"))
      {

      typename TubeType::Pointer tub = TubeType::New();
      tubespacing = (*TubeIterator)->GetSpacing();

      tub->SetId(((TubeType *)((*TubeIterator).GetPointer()))->GetId());
      tub->SetRoot(((TubeType *)((*TubeIterator).GetPointer()))->GetRoot());
      tub->SetArtery(((TubeType *)((*TubeIterator).GetPointer()))
        ->GetArtery());

      tub->GetProperty()->SetColor((*TubeIterator)->GetProperty()
        ->GetColor());

      TubePointListType tubeList =
        ((TubeType *)((*TubeIterator).GetPointer()))->GetPoints();
      typename TubePointListType::const_iterator TubePointIterator =
        tubeList.begin();

      while(TubePointIterator != tubeList.end())
        {
        point = (*TubePointIterator).GetPosition();
        for(unsigned int i=0;i<TDimension;i++)
          {
          point[i] *= tubespacing[i]*groupspacing[i]*m_Scale;
          }

        if(m_TransformAsGroup)
          {
          point = m_TransformAsGroup->GetObjectToParentTransform()->TransformPoint(point);
          }
        else
          {
          point = m_Transform->TransformPoint(point);
          }

        for(unsigned int i=0;i<TDimension;i++)
          {
          point[i] /= tubespacing[i]*groupspacing[i];
          }

        /** Crop the tube net to fit the image*/
        bool IsInside = true;
        if(m_Crop)
          {
          for(unsigned int i=0;i<TDimension;i++)
            {
            if( ( point[i] > m_CropSize[i] + m_NarrowBandSize )
                || ( point[i] < 0 ) )
                // no negative numbers, the transformation should
                // take care of the narrrowband
              {
              IsInside = false;
              break;
              }
            }
          }

        if(IsInside)
          {
          VesselTubeSpatialObjectPoint<TDimension> pnt;
          pnt.SetPosition(point);

          if(m_TransformAsGroup)
            {
            FirstNormal = m_TransformAsGroup->GetObjectToParentTransform()
              ->TransformCovariantVector(
                (*TubePointIterator).GetNormal1() );
            SecondNormal  = m_TransformAsGroup->GetObjectToParentTransform()
              ->TransformCovariantVector(
                (*TubePointIterator).GetNormal2() );
            pnt.SetNormal1(FirstNormal);
            pnt.SetNormal2(SecondNormal);
            }
          else
            {
            FirstNormal = m_Transform->TransformCovariantVector(
              (*TubePointIterator).GetNormal1() );
            SecondNormal  = m_Transform->TransformCovariantVector(
              (*TubePointIterator).GetNormal2() );
            pnt.SetNormal1(FirstNormal);
            pnt.SetNormal2(SecondNormal);
            }

          pnt.SetRadius((*TubePointIterator).GetRadius()*m_Scale);

          if(m_Medialness == -1)
            {
            pnt.SetMedialness((*TubePointIterator).GetMedialness());
            }
          else
            {
            pnt.SetMedialness(m_Medialness);
            }

          if(m_Ridgeness == -1)
            {
            pnt.SetRidgeness((*TubePointIterator).GetRidgeness());
            }
          else
            {
            pnt.SetRidgeness(m_Ridgeness);
            }

          pnt.SetBranchness((*TubePointIterator).GetBranchness());
          tub->GetPoints().push_back(pnt);
          }
        TubePointIterator++;
        }

      tub->GetIndexToObjectTransform()->SetScaleComponent(
        ((*TubeIterator)->GetIndexToObjectTransform()
          ->GetScaleComponent()));
      tub->RemoveDuplicatePoints();
      tub->ComputeTangentAndNormals();
      m_Output->AddSpatialObject(tub);
      }
    }
  delete TubeList;
  delete [] groupspacing;
}


template <class TTransformType, unsigned int TDimension>
void
TubeToTubeTransformFilter<TTransformType,TDimension>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Transformation: " << m_Transform << std::endl;
}

} // End namespace itk

#endif // End !defined(__itkTubeToTubeTransformFilter_txx)
