/*=========================================================================

Library:   TubeTK/VTree3D

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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

#ifndef __itktubeTubeExtractor_hxx
#define __itktubeTubeExtractor_hxx


#include <itktubeLimitedMinimumMaximumImageFilter.h>

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkImageDuplicator.h>

namespace itk
{

namespace tube
{

/**
 * Constructor */
template <class TInputImage>
TubeExtractor<TInputImage>::TubeExtractor(void)
{
  m_RidgeExtractor = RidgeExtractorType::New();
  m_RadiusExtractor = RadiusExtractorType::New();
  m_RidgeExtractor->SetRadiusExtractor(m_RadiusExtractor);

  m_OptimizeRadius = true;

  m_IdleCallBack = nullptr;
  m_StatusCallBack = nullptr;
  m_NewTubeCallBack = nullptr;
  m_AbortProcess = nullptr;

  m_TubeGroup = TubeGroupType::New();

  m_TubeColor.set_size(4);
  m_TubeColor[0] = 1.0f;
  m_TubeColor[1] = 0.0f;
  m_TubeColor[2] = 0.0f;
  m_TubeColor[3] = 1.0f;

  m_SeedMask = nullptr;
  m_SeedRadiusMask = nullptr;
  m_SeedMaskMaximumNumberOfPoints = 0;
  m_SeedMaskStride = 1;
  m_UseSeedMaskAsProbabilities = false;
  m_SeedExtractionMinimumSuccessRatio = 0;
  m_SeedExtractionMinimumProbability = -9999999;

  m_SeedsInObjectSpaceList.clear();
  m_SeedRadiiInObjectSpaceList.clear();
}

/**
 * Destructor */
template <class TInputImage>
TubeExtractor<TInputImage>::~TubeExtractor(void)
{
  m_SeedsInObjectSpaceList.clear();
  m_SeedRadiiInObjectSpaceList.clear();
}

/**
 * Set the input image */
template <class TInputImage>
void
TubeExtractor<TInputImage>::SetInputImage(ImageType * inputImage)
{
  // this->m_RidgeExtractor = RidgeExtractor<ImageType>::New();
  this->m_RidgeExtractor->SetInputImage(inputImage);

  // this->m_RadiusExtractor = RadiusExtractor3<ImageType>::New();
  this->m_RadiusExtractor->SetInputImage(inputImage);
}

template <class TInputImage>
const TInputImage *
TubeExtractor<TInputImage>::GetInputImage(void) const
{
  return this->m_RidgeExtractor->GetInputImage();
}

/**
 * Optionally set a different image for radius estimation */
template <class TInputImage>
void
TubeExtractor<TInputImage>::SetRadiusInputImage(ImageType * inputImage)
{
  this->m_RadiusExtractor->SetInputImage(inputImage);

  this->m_RidgeExtractor->SetRadiusExtractor(this->m_RadiusExtractor);
}

template <class TInputImage>
const TInputImage *
TubeExtractor<TInputImage>::GetRadiusInputImage(void) const
{
  return this->m_RadiusExtractor->GetInputImage();
}


/**
 * Set the tube mask image */
template <class TInputImage>
void
TubeExtractor<TInputImage>::SetTubeMaskImage(typename TubeExtractor<TInputImage>::TubeMaskImageType * mask)
{
  m_RidgeExtractor->SetTubeMaskImage(mask);
}

template <class TInputImage>
typename TubeExtractor<TInputImage>::TubeMaskImageType *
TubeExtractor<TInputImage>::GetTubeMaskImage(void)
{
  return this->m_RidgeExtractor->GetTubeMaskImage();
}


/**
 * Set Data Min value */
template <class TInputImage>
void
TubeExtractor<TInputImage>::SetDataMin(double dataMin)
{
  if (this->m_RidgeExtractor.IsNull() || this->m_RadiusExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  this->m_RidgeExtractor->SetDataMin(dataMin);
  this->m_RadiusExtractor->SetDataMin(dataMin);
}

/**
 * Get Data Min value */
template <class TInputImage>
double
TubeExtractor<TInputImage>::GetDataMin(void)
{
  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  return this->m_RidgeExtractor->GetDataMin();
}


template <class TInputImage>
void
TubeExtractor<TInputImage>::SetDataMinMaxLimits(double limitMin, double limitMax)
{
  typedef LimitedMinimumMaximumImageFilter<ImageType> MinMaxFilterType;
  typename MinMaxFilterType::Pointer                  minMaxFilter = MinMaxFilterType::New();
  minMaxFilter->SetInput(this->GetInputImage());
  minMaxFilter->SetMinimumLimit(limitMin);
  minMaxFilter->SetMaximumLimit(limitMax);
  minMaxFilter->Update();
  this->SetDataMin(minMaxFilter->GetMinimum());
  this->SetDataMax(minMaxFilter->GetMaximum());
}

/**
 * Set Data Max value */
template <class TInputImage>
void
TubeExtractor<TInputImage>::SetDataMax(double dataMax)
{
  if (this->m_RidgeExtractor.IsNull() || this->m_RadiusExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  this->m_RidgeExtractor->SetDataMax(dataMax);
  this->m_RadiusExtractor->SetDataMax(dataMax);
}


/**
 * Get Data Max value */
template <class TInputImage>
double
TubeExtractor<TInputImage>::GetDataMax(void)
{
  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  return this->m_RidgeExtractor->GetDataMax();
}

/**
 * Set border value */
template <class TInputImage>
void
TubeExtractor<TInputImage>::SetBorderInIndexSpace(int border)
{
  typename ImageType::IndexType minIndx = this->GetInputImage()->GetLargestPossibleRegion().GetIndex();
  typename ImageType::SizeType  size = this->GetInputImage()->GetLargestPossibleRegion().GetSize();
  typename ImageType::IndexType maxIndx;
  for (unsigned int i = 0; i < ImageDimension; ++i)
  {
    maxIndx[i] = minIndx[i] + size[i] - 1 - border;
    minIndx[i] += border;
  }
  this->SetExtractBoundMinInIndexSpace(minIndx);
  this->SetExtractBoundMaxInIndexSpace(maxIndx);
}


template <class TInputImage>
void
TubeExtractor<TInputImage>::SetExtractBoundMinInIndexSpace(const typename TInputImage::IndexType & dataMin)
{
  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  this->m_RidgeExtractor->SetExtractBoundMinInIndexSpace(dataMin);
}

/**
 * Get Data Min value */
template <class TInputImage>
typename TInputImage::IndexType
TubeExtractor<TInputImage>::GetExtractBoundMinInIndexSpace(void) const
{
  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  return this->m_RidgeExtractor->GetExtractBoundMinInIndexSpace();
}

/**
 * Set Data Max value */
template <class TInputImage>
void
TubeExtractor<TInputImage>::SetExtractBoundMaxInIndexSpace(const typename TInputImage::IndexType & boundMax)
{
  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  this->m_RidgeExtractor->SetExtractBoundMaxInIndexSpace(boundMax);
}


/**
 * Get Data Max value */
template <class TInputImage>
typename TInputImage::IndexType
TubeExtractor<TInputImage>::GetExtractBoundMaxInIndexSpace(void) const
{
  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  return this->m_RidgeExtractor->GetExtractBoundMaxInIndexSpace();
}

/**
 * Set Radius */
template <class TInputImage>
void
TubeExtractor<TInputImage>::SetRadiusInObjectSpace(double radius)
{
  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  this->m_RidgeExtractor->SetScale(radius);
  this->m_RadiusExtractor->SetRadiusStart(radius);
}

/**
 * Get Radius */
template <class TInputImage>
double
TubeExtractor<TInputImage>::GetRadiusInObjectSpace(void)
{
  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  return this->m_RidgeExtractor->GetScale();
}

/**
 * Get the ridge extractor */
template <class TInputImage>
RidgeExtractor<TInputImage> *
TubeExtractor<TInputImage>::GetRidgeExtractor(void)
{
  return this->m_RidgeExtractor.GetPointer();
}

/**
 * Get the radius extractor */
template <class TInputImage>
RadiusExtractor3<TInputImage> *
TubeExtractor<TInputImage>::GetRadiusExtractor(void)
{
  return this->m_RadiusExtractor.GetPointer();
}

/**
 * Extract the tube given the position of the first point
 * and the tube ID */
template <class TInputImage>
bool
TubeExtractor<TInputImage>::FindLocalTubeInObjectSpace(PointType & x)
{
  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  return this->m_RidgeExtractor->LocalRidge(x);
}

/**
 * Extract the tube given the position of the first point
 * and the tube ID */
template <class TInputImage>
typename TubeExtractor<TInputImage>::TubeType *
TubeExtractor<TInputImage>::ExtractTubeInObjectSpace(const PointType & x, unsigned int tubeID, bool verbose)
{
  if (verbose)
  {
    std::cout << "TubeExtractor: ExtracTubeInObjectSpace: Start" << std::endl;
  }

  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  IndexType xi;
  if (!this->m_RidgeExtractor->GetTubeMaskImage()->TransformPhysicalPointToIndex(x, xi))
  {
    if (verbose)
    {
      std::cout << "Point maps to outside of image. Aborting." << std::endl;
      return nullptr;
    }
  }

  if (verbose)
  {
    std::cout << "Physical point = " << x << std::endl;
    std::cout << "Index point = " << xi << std::endl;
    std::cout << "Mask value = " << this->m_RidgeExtractor->GetTubeMaskImage()->GetPixel(xi) << std::endl;
  }

  if (this->m_RidgeExtractor->GetTubeMaskImage()->GetPixel(xi) != 0)
  {
    if (verbose || this->GetDebug())
    {
      std::cout << "Initial pixel on prior tube." << std::endl;
      std::cout << "  x = " << x << std::endl;
      std::cout << "  xi = " << xi << std::endl;
    }
    return nullptr;
  }
  else if (verbose)
  {
    std::cout << "No overlapping tube" << std::endl;
  }

  typename TubeType::Pointer tube = this->m_RidgeExtractor->ExtractRidge(x, tubeID, verbose);

  if (tube.IsNull())
  {
    if (verbose || this->GetDebug())
    {
      std::cout << "m_RidgeExtractor->Extract() fails!" << std::endl;
      std::cout << "  x = " << x << std::endl;
    }
    return nullptr;
  }

  if (this->m_AbortProcess != NULL)
  {
    if (this->m_AbortProcess())
    {
      if (this->m_StatusCallBack)
      {
        this->m_StatusCallBack("Extract: Ridge", "Aborted", 0);
      }
      return nullptr;
    }
  }

  if (m_OptimizeRadius)
  {
    if (!this->m_RadiusExtractor->ExtractRadii(tube, verbose))
    {
      return nullptr;
    }
  }
  else
  {
    if (m_SeedRadiusMask)
    {
      double                                        defaultR = this->m_RadiusExtractor->GetRadiusStart();
      typename std::vector<TubePointType>::iterator pntIter;
      pntIter = tube->GetPoints().begin();
      typename std::vector<TubePointType>::iterator pntIterEnd;
      pntIterEnd = tube->GetPoints().end();
      while (pntIter != pntIterEnd)
      {
        PointType pnt = pntIter->GetPositionInObjectSpace();
        IndexType indx;
        if (m_SeedRadiusMask->TransformPhysicalPointToIndex(pnt, indx))
        {
          double r = m_SeedRadiusMask->GetPixel(indx);
          if (r != 0)
          {
            pntIter->SetRadiusInObjectSpace(r);
          }
          else
          {
            pntIter->SetRadiusInObjectSpace(defaultR);
          }
        }
        ++pntIter;
      }
    }
  }


  if (this->m_NewTubeCallBack != NULL)
  {
    this->m_NewTubeCallBack(tube);
  }

  if (this->m_StatusCallBack)
  {
    char s[80];
    std::snprintf(s, 80, "%zd points", tube->GetPoints().size());
    this->m_StatusCallBack("Extract: Ridge", s, 0);
  }

  if (verbose)
  {
    std::cout << "Adding tube to group." << std::endl;
  }
  this->AddTube(tube);

  tube->Register();
  return tube;
}

template <class TInputImage>
void
TubeExtractor<TInputImage>::SetSeedsInIndexSpaceList(const ContinuousIndexListType & iList)
{
  m_SeedsInObjectSpaceList.clear();
  m_SeedRadiiInObjectSpaceList.clear();
  double    radiusInObjectSpace = this->m_RadiusExtractor->GetRadiusStart();
  PointType pnt;
  for (size_t seedNum = 0; seedNum < iList.size(); ++seedNum)
  {
    this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(iList[seedNum], pnt);
    m_SeedsInObjectSpaceList.push_back(pnt);
    m_SeedRadiiInObjectSpaceList.push_back(radiusInObjectSpace);
  }
}

template <class TInputImage>
void
TubeExtractor<TInputImage>::SetSeedsInObjectSpaceList(const PointListType & oList)
{
  m_SeedsInObjectSpaceList.clear();
  m_SeedRadiiInObjectSpaceList.clear();
  double radiusInObjectSpace = this->m_RadiusExtractor->GetRadiusStart();
  for (size_t seedNum = 0; seedNum < oList.size(); ++seedNum)
  {
    m_SeedsInObjectSpaceList.push_back(oList[seedNum]);
    m_SeedRadiiInObjectSpaceList.push_back(radiusInObjectSpace);
  }
}

template <class TInputImage>
void
TubeExtractor<TInputImage>::SetSeedRadiiInObjectSpaceList(const RadiusListType & rList)
{
  m_SeedRadiiInObjectSpaceList.clear();
  for (size_t seedNum = 0; seedNum < rList.size(); ++seedNum)
  {
    m_SeedRadiiInObjectSpaceList.push_back(rList[seedNum]);
  }
}

template <class TInputImage>
void
TubeExtractor<TInputImage>::ProcessSeeds(bool verbose)
{
  this->GetRidgeExtractor()->ResetFailureCodeCounts();
  double defaultR = this->GetRadiusInObjectSpace();

  if (this->m_SeedMask.IsNotNull())
  {
    if (m_UseSeedMaskAsProbabilities)
    {
      typedef itk::ImageDuplicator<TubeMaskImageType> DuplicatorType;
      typename DuplicatorType::Pointer                duplicator = DuplicatorType::New();
      duplicator->SetInputImage(m_SeedMask);
      duplicator->Update();
      typename TubeMaskImageType::Pointer tmpSeedMask = duplicator->GetOutput();

      typedef itk::MinimumMaximumImageCalculator<TubeMaskImageType> MinMaxCalcType;
      typename MinMaxCalcType::Pointer                              maxCalc = MinMaxCalcType::New();
      unsigned int                                                  count = 1;
      double                                                        successRatio = 1;
      double                                                        maxValue = m_SeedExtractionMinimumProbability;
      while ((m_SeedMaskMaximumNumberOfPoints == 0 || count < m_SeedMaskMaximumNumberOfPoints) &&
             successRatio >= m_SeedExtractionMinimumSuccessRatio && maxValue >= m_SeedExtractionMinimumProbability)
      {
        std::cout << "Count = " << count << std::endl;
        maxCalc->SetImage(tmpSeedMask);
        maxCalc->ComputeMaximum();
        maxValue = maxCalc->GetMaximum();
        typename ImageType::IndexType maxIndx = maxCalc->GetIndexOfMaximum();
        if (maxValue >= m_SeedExtractionMinimumProbability)
        {
          if (this->m_SeedRadiusMask)
          {
            double r = m_SeedRadiusMask->GetPixel(maxIndx);
            if (r <= 0)
            {
              r = defaultR;
            }
            this->SetRadiusInObjectSpace(r);
          }
          PointType pnt;
          this->GetInputImage()->TransformIndexToPhysicalPoint(maxIndx, pnt);
          typename TubeType::Pointer xTube = this->ExtractTubeInObjectSpace(pnt, count, verbose);
          if (!xTube.IsNull())
          {
            m_RidgeExtractor->DeleteTube(xTube, tmpSeedMask.GetPointer());
            successRatio = (successRatio * 9 + 1) / 10;
            std::cout << "   Ridge size = " << xTube->GetNumberOfPoints() << std::endl;
          }
          else
          {
            typedef itk::NeighborhoodIterator<TubeMaskImageType> NeighborIterType;
            typename NeighborIterType::RadiusType                radius;
            radius.Fill(5);
            NeighborIterType iter(radius, tmpSeedMask, tmpSeedMask->GetLargestPossibleRegion());
            iter.SetLocation(maxIndx);
            bool inside;
            for (unsigned int i = 0; i < iter.Size(); ++i)
            {
              iter.SetPixel(i, 0, inside);
            }
            successRatio = (successRatio * 9 + 0) / 10;
            std::cout << "   Ridge not found" << std::endl;
          }
        }
        ++count;
      }
    }
    else
    {
      ImageRegionConstIteratorWithIndex<TubeMaskImageType> iter(this->GetSeedMask(),
                                                                this->GetSeedMask()->GetLargestPossibleRegion());
      ImageRegionConstIterator<ImageType>                  iterR;
      if (this->m_SeedRadiusMask)
      {
        iterR = ImageRegionConstIterator<ImageType>(this->m_SeedRadiusMask,
                                                    this->m_SeedRadiusMask->GetLargestPossibleRegion());
      }

      int    count = 0;
      double radiusInObjectSpace = this->m_RadiusExtractor->GetRadiusStart();
      while (!iter.IsAtEnd())
      {
        if (iter.Get())
        {
          if (++count == this->m_SeedMaskStride)
          {
            count = 0;
            PointType pnt;
            this->GetInputImage()->TransformIndexToPhysicalPoint(iter.GetIndex(), pnt);
            m_SeedsInObjectSpaceList.push_back(pnt);
            if (this->m_SeedRadiusMask.IsNotNull())
            {
              m_SeedRadiiInObjectSpaceList.push_back(iterR.Get());
              ++iterR;
            }
            else
            {
              m_SeedRadiiInObjectSpaceList.push_back(radiusInObjectSpace);
            }
          }
        }
        ++iter;
      }
    }
  }

  if (m_SeedsInObjectSpaceList.size() > 0)
  {
    typename std::vector<PointType>::iterator seedIter = this->m_SeedsInObjectSpaceList.begin();
    typename std::vector<double>::iterator    seedRadiusIter = this->m_SeedRadiiInObjectSpaceList.begin();

    bool useRadiiList = false;
    if (m_SeedRadiiInObjectSpaceList.size() == m_SeedsInObjectSpaceList.size())
    {
      useRadiiList = true;
    }

    unsigned int count = 1;
    unsigned int maxCount = m_SeedsInObjectSpaceList.size();
    bool         foundOneTube = false;
    while (seedIter != this->m_SeedsInObjectSpaceList.end())
    {
      PointType x = *seedIter;

      std::cout << "Extracting from index point " << x << " (" << (count / (double)maxCount) * 100 << "%)" << std::endl;

      if (useRadiiList)
      {
        this->SetRadiusInObjectSpace(*seedRadiusIter);
        ++seedRadiusIter;
      }

      typename TubeType::Pointer xTube = this->ExtractTubeInObjectSpace(x, count, verbose);
      if (!xTube.IsNull())
      {
        foundOneTube = true;
        std::cout << "   Ridge size = " << xTube->GetNumberOfPoints() << std::endl;
      }
      else
      {
        std::cout << "   Ridge not found" << std::endl;
      }

      ++seedIter;
      ++count;
    }
    if (!foundOneTube)
    {
      std::cout << "*** No Ridges found! ***" << std::endl;
      return;
    }
  }

  std::cout << "Ridge termination code counts:" << std::endl;
  for (unsigned int code = 0; code < this->GetRidgeExtractor()->GetNumberOfFailureCodes(); ++code)
  {
    std::cout << "   " << this->m_RidgeExtractor->GetFailureCodeName(typename RidgeExtractorType::FailureCodeEnum(code))
              << " : "
              << this->m_RidgeExtractor->GetFailureCodeCount(typename RidgeExtractorType::FailureCodeEnum(code))
              << std::endl;
  }
}

/**
 * Get list of extracted tubes */
template <class TInputImage>
typename TubeExtractor<TInputImage>::TubeGroupType *
TubeExtractor<TInputImage>::GetTubeGroup(void)
{
  return m_TubeGroup;
}

/**
 * Set list of extracted tubes */
template <class TInputImage>
void
TubeExtractor<TInputImage>::SetTubeGroup(TubeGroupType * tubes)
{
  m_TubeGroup = tubes;
  typename TubeGroupType::ChildrenListType *         cList = tubes->GetChildren(9999);
  typename TubeGroupType::ChildrenListType::iterator iter = cList->begin();
  while (iter != cList->end())
  {
    this->AddTube(static_cast<TubeType *>(iter->GetPointer()));
    ++iter;
  }
}

/**
 * Smooth a tube */
template <class TInputImage>
void
TubeExtractor<TInputImage>::SmoothTube(TubeType * tube, int h)
{
  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  ::tube::TubeMathFilters<ImageDimension> filter;
  filter.SetInputTube(tube);
  filter.SmoothTube(h);
}

/**
 * Add a tube */
template <class TInputImage>
bool
TubeExtractor<TInputImage>::AddTube(TubeType * tube)
{
  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  bool result = this->m_RidgeExtractor->AddTube(tube);
  if (result)
  {
    m_TubeGroup->AddChild(tube);
  }

  return result;
}

/**
 * Delete a tube */
template <class TInputImage>
bool
TubeExtractor<TInputImage>::DeleteTube(TubeType * tube)
{
  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  bool result = this->m_RidgeExtractor->DeleteTube(tube);
  if (result)
  {
    m_TubeGroup->RemoveChild(tube);
  }

  return result;
}

/**
 * Set the tube color */
template <class TInputImage>
void
TubeExtractor<TInputImage>::SetTubeColor(const vnl_vector<double> & color)
{
  int nc = color.size();
  if (nc > 4)
  {
    nc = 4;
  }
  else if (nc < 4)
  {
    this->m_TubeColor[3] = 1.0;
  }
  for (int i = 0; i < nc; i++)
  {
    this->m_TubeColor[i] = color[i];
  }
}

template <class TInputImage>
vnl_vector<double> &
TubeExtractor<TInputImage>::GetTubeColor(void)
{
  return m_TubeColor;
}

/**
 * Set the idle call back */
template <class TInputImage>
void
TubeExtractor<TInputImage>::IdleCallBack(bool (*idleCallBack)())
{
  this->m_IdleCallBack = idleCallBack;
}

/**
 * Set the status callback  */
template <class TInputImage>
void
TubeExtractor<TInputImage>::StatusCallBack(void (*statusCallBack)(const char *, const char *, int))
{
  if (this->m_RidgeExtractor.IsNull())
  {
    throw("Input data must be set first in TubeExtractor");
  }

  this->m_StatusCallBack = statusCallBack;
  this->m_RidgeExtractor->StatusCallBack(statusCallBack);
  this->m_RadiusExtractor->StatusCallBack(statusCallBack);
}

/**
 * Set the status callback  */
template <class TInputImage>
void
TubeExtractor<TInputImage>::NewTubeCallBack(void (*newTubeCallBack)(TubeType *))
{
  this->m_NewTubeCallBack = newTubeCallBack;
}

/**
 * Abort the process  */
template <class TInputImage>
void
TubeExtractor<TInputImage>::AbortProcess(bool (*abortProcess)())
{
  this->m_AbortProcess = abortProcess;
}

/**
 * PrintSelf */
template <class TInputImage>
void
TubeExtractor<TInputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "RidgeExtractor = " << this->m_RidgeExtractor << std::endl;
  os << indent << "RadiusExtractor = " << this->m_RadiusExtractor << std::endl;

  os << indent << "TubeGroup = " << this->m_TubeGroup << std::endl;
  os << indent << "SeedsInObjectSpaceList.size = " << this->m_SeedsInObjectSpaceList.size() << std::endl;
  os << indent << "SeedRadiiInObjectSpaceList.size = " << this->m_SeedRadiiInObjectSpaceList.size() << std::endl;
  os << indent << "SeedMask = " << this->m_SeedMask << std::endl;
  os << indent << "SeedRadiusMask = " << this->m_SeedRadiusMask << std::endl;
  os << indent << "SeedMaskStride = " << this->m_SeedMaskStride << std::endl;

  os << indent << "TubeColor.r = " << this->m_TubeColor[0] << std::endl;
  os << indent << "TubeColor.g = " << this->m_TubeColor[1] << std::endl;
  os << indent << "TubeColor.b = " << this->m_TubeColor[2] << std::endl;
  os << indent << "TubeColor.a = " << this->m_TubeColor[3] << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeTubeExtractor_hxx )
