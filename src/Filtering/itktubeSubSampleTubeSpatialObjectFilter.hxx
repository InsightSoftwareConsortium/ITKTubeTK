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

#ifndef __itktubeSubSampleTubeSpatialObjectFilter_hxx
#define __itktubeSubSampleTubeSpatialObjectFilter_hxx


namespace itk
{

namespace tube
{

template <unsigned int ObjectDimension>
SubSampleTubeSpatialObjectFilter<ObjectDimension>::SubSampleTubeSpatialObjectFilter(void)
  : m_Sampling(1)
{}

template <unsigned int ObjectDimension>
SubSampleTubeSpatialObjectFilter<ObjectDimension>::~SubSampleTubeSpatialObjectFilter(void)
{}


template <unsigned int ObjectDimension>
void
SubSampleTubeSpatialObjectFilter<ObjectDimension>::GenerateData(void)
{
  const TubeSpatialObjectType * input = dynamic_cast<const TubeSpatialObjectType *>(this->GetInput());
  if (input == nullptr)
  {
    std::cerr << "Error: tube passed to SubSampleTubes is not a tube." << std::endl;
    return;
  }
  typename TubeSpatialObjectType::Pointer output = dynamic_cast<TubeSpatialObjectType *>(this->GetOutput());

  typedef typename TubeSpatialObjectType::TubePointListType TubePointListType;
  const TubePointListType &                                 inputPoints = input->GetPoints();
  TubePointListType &                                       outputPoints = output->GetPoints();

  const unsigned int numberOfInputPoints = inputPoints.size();
  unsigned int       numberOfOutputPoints;
  if (this->m_Sampling == 1)
  {
    numberOfOutputPoints = numberOfInputPoints / this->m_Sampling + 0;
  }
  else if (numberOfInputPoints % this->m_Sampling == 0)
  {
    numberOfOutputPoints = numberOfInputPoints / this->m_Sampling + 1;
  }
  else
  {
    numberOfOutputPoints = numberOfInputPoints / this->m_Sampling + 2;
  }
  outputPoints.resize(numberOfOutputPoints);
  for (unsigned int inputIndex = 0, outputIndex = 0; outputIndex < numberOfOutputPoints - 1;
       ++outputIndex, inputIndex += this->m_Sampling)
  {
    outputPoints[outputIndex] = inputPoints[inputIndex];
  }
  outputPoints[numberOfOutputPoints - 1] = inputPoints[numberOfInputPoints - 1];

  output->RemoveDuplicatePointsInObjectSpace();

  this->GraftOutput(output);
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeSubSampleTubeSpatialObjectFilter_hxx )
