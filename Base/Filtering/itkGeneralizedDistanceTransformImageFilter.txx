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

#ifndef __itkGeneralizedDistanceTransformImageFilter_txx
#define __itkGeneralizedDistanceTransformImageFilter_txx

#include <limits>

#include "itkGeneralizedDistanceTransformImageFilter.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace
{
// The intersection method needs to divide by a value i*s^2 for an integer i
// with 0 < i < maximal image extent, and a spacing s.  We calculate a table
// of reciprocal values and multiply by those instead to improve speed.
//
// It is a little less precise, but in the context of parabola intersection it
// makes hardly a difference.
template <class TSpacingType>
void updateDivisionTable(
  std::vector< TSpacingType > & table, TSpacingType spacing, size_t size)
{
  for (size_t i = 1; i < size; ++i)
    {
    table[i] = TSpacingType(1)/(spacing*spacing*i);
    }
}
}

namespace itk
{

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
void
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::SetCreateVoronoiMap(bool b)
{
  m_CreateVoronoiMap = b;
  if( m_CreateVoronoiMap )
    {
    this->SetNumberOfRequiredInputs(2);
    this->SetNumberOfRequiredOutputs(2);
    }
  else
    {
    this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs(1);
    }
  this->Modified();
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
void
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::CreateVoronoiMapOn()
{
  this->SetCreateVoronoiMap(true);
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
void
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::CreateVoronoiMapOff()
{
  this->SetCreateVoronoiMap(false);
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
double
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::GetMaximalSquaredDistance() const
{
  return m_MaximalSquaredDistance;
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
void
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::SetInput1(const FunctionImageType *functionImage)
{
  // Process object is not const-correct so the const casting is required.
  this->SetNthInput(0, const_cast<FunctionImageType *>(functionImage));
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
void
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::SetInput2(const LabelImageType *labelImage)
{
  this->SetNthInput(1, const_cast<LabelImageType *>(labelImage));
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
typename
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::DistanceImageType*
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::GetDistance(void)
{
  return  dynamic_cast<DistanceImageType *>(this->ProcessObject::GetOutput(0));
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
typename
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::LabelImageType*
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::GetVoronoiMap()
{
  assert(this->m_CreateVoronoiMap);
  return  dynamic_cast<LabelImageType *>(this->ProcessObject::GetOutput(1));
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::GeneralizedDistanceTransformImageFilter()
  : m_MaximalSquaredDistance(
      std::numeric_limits<typename TDistanceImage::PixelType>::max()/2),
    m_CacheLineSize(128)
{
  // First check the constraints on the types
  // All numeric types have to be signed. It can be argued that
  // this is overly restrictive for DistancePixelType when you only use a
  // non-negative input function, but internally a negative value of type
  // DistancePixelType is used during sampling.
  assert(std::numeric_limits<AbscissaIndexType>::is_signed);
  assert(std::numeric_limits<DistancePixelType>::is_signed);
  assert(std::numeric_limits<SpacingType>::is_signed);

  // The interval covered by SpacingType must contain those covered by
  // AbscissaIndexType and ApexHeightType. We only check for the maximum.
  assert(std::numeric_limits<SpacingType>::max() >=
         std::numeric_limits<AbscissaIndexType>::max());
  assert(std::numeric_limits<SpacingType>::max() >=
         std::numeric_limits<DistancePixelType>::max());

  // The filter defaults to being fast.
  this->CreateVoronoiMapOff();
  this->UseImageSpacingOff();

  // Create output pointers. Data is allocated later, so it does not hurt to
  // create a voronoi map right here, even if it will not be filled with data by
  // default.
  DistanceImagePointer distance = DistanceImageType::New();
  this->SetNthOutput(0, distance.GetPointer());

  LabelImagePointer voronoiMap = LabelImageType::New();
  this->SetNthOutput(1, voronoiMap.GetPointer());
}

/**
 *  Print Self
 */
template < class TFunctionImage,class TDistanceImage, class TLabelImage >
void
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Generalized Distance Transform: " << std::endl;
  os << indent << "MaximalSquaredDistance: "
    << m_MaximalSquaredDistance << std::endl;
  os << indent << "UseSpacing: " << this->m_UseImageSpacing << std::endl;
  os << indent << "CreateVoronoiMap: " << this->m_CreateVoronoiMap << std::endl;
}

template< class TFunctionImage, class TDistanceImage, class TLabelImage >
void
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::EnlargeOutputRequestedRegion(DataObject *)
{
  this->GetDistance()->SetRequestedRegion(
    this->GetDistance()->GetLargestPossibleRegion());
  if (this->m_CreateVoronoiMap)
    {
    this->GetVoronoiMap()->SetRequestedRegion(
      this->GetVoronoiMap()->GetLargestPossibleRegion());
    }
}

/**
 * Allocate and initialize output images. Helper function for GenerateData()
 */
template < class TFunctionImage,class TDistanceImage, class TLabelImage >
void
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::PrepareData()
{
  // Copy the function image into the distance image
  FunctionImageConstPointer functionImage  =
    dynamic_cast<FunctionImageType *>(ProcessObject::GetInput(0));

  std::cout << "spacing = " << functionImage->GetSpacing() << std::endl;
  DistanceImagePointer distance = this->GetDistance();
  distance->SetRegions( functionImage->GetLargestPossibleRegion() );
  distance->CopyInformation( functionImage );
  std::cout << "spacing = " << distance->GetSpacing() << std::endl;
  distance->Allocate();

  ImageRegionConstIterator<FunctionImageType>
    functionIt(functionImage, distance->GetRequestedRegion());
  ImageRegionIterator<DistanceImageType>
    distanceIt(distance, distance->GetRequestedRegion());

  functionIt.GoToBegin();
  distanceIt.GoToBegin();
  while(!distanceIt.IsAtEnd())
    {
    // Use functionIt.Get() in order to account for image adaptors. distanceIt
    // is our own, so we can write to it using distanceIt.Value().
    distanceIt.Value() =
      static_cast<typename DistanceImageType::PixelType>(functionIt.Get());

    ++functionIt;
    ++distanceIt;
    }

  if (this->m_CreateVoronoiMap)
    {
    // Copy the label image into the voronoi map
    LabelImagePointer labelImage  =
      dynamic_cast<LabelImageType *>(ProcessObject::GetInput(1));

    LabelImagePointer voronoiMap = this->GetVoronoiMap();
    voronoiMap->SetRegions( labelImage->GetLargestPossibleRegion() );
    // voronoiMap->SetBufferedRegion(voronoiMap->GetRequestedRegion());
    voronoiMap->CopyInformation( labelImage );
    voronoiMap->Allocate();

    ImageRegionConstIterator<LabelImageType>
      labelIt(labelImage, distance->GetRequestedRegion());
    ImageRegionIterator<LabelImageType>
      voronoiIt(voronoiMap, distance->GetRequestedRegion());

    labelIt.GoToBegin();
    voronoiIt.GoToBegin();
    while(!voronoiIt.IsAtEnd())
      {
      voronoiIt.Value() = labelIt.Get();
      ++labelIt;
      ++voronoiIt;
      }
    }
}

/**
 *  Compute Distance and Voronoi maps
 *  \todo Support progress methods/callbacks.
 */
template < class TFunctionImage,class TDistanceImage, class TLabelImage >
void
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::GenerateData()
{
  this->PrepareData();

  // We need the size and probably the spacing of the images.
  DistanceImagePointer distance = this->GetDistance();
  typename DistanceImageType::SpacingType spacing = distance->GetSpacing();
  typename DistanceImageType::SizeType size =
    distance->GetRequestedRegion().GetSize();

  // The distance image has been initialized to contain the function values
  // f(x) at x = (x1 x2 ... xN).
  // It is transformed into the lower envelope of spherical paraboloids rooted
  // at (x1 x2 ... xN f(x)) by iteration over the dimensions of the image.
  // Information on the region covered by a paraboloid is provided optionally
  // by copying the label at x.
  //
  // The iterations visit each scanline in each dimension.
  //
  // To make better use of L1 cache, a few scanlines are read and written in
  // parallel for dimensions 1 and upward. See below for details.
  typedef itk::ImageLinearIteratorWithIndex<DistanceImageType> DIt;
  DIt distanceIt(distance, distance->GetRequestedRegion());

  typedef itk::ImageLinearIteratorWithIndex<LabelImageType> LIt;
  LIt voronoiMapIt;
  if (this->m_CreateVoronoiMap)
    {
    LabelImagePointer voronoiMap =
      dynamic_cast<LabelImageType *>(this->ProcessObject::GetOutput(1));
    voronoiMapIt = LIt(voronoiMap, voronoiMap->GetRequestedRegion());
    }

  // Compute the maximal extent in number of pixels, rounded up to fill full
  // cache lines.
  // This is used to set up various vectors, buffers, and a division table
  typedef typename DistanceImageType::SizeValueType DistSizeValueType;
  DistSizeValueType maxSize = 0;

  for (unsigned int d = 0; d < FunctionImageType::ImageDimension; ++d)
    {
    maxSize = std::max(maxSize, size[d]);
    }

  maxSize = (DistSizeValueType)(
    std::ceil((double)maxSize / m_CacheLineSize) * m_CacheLineSize);

  // With the division table, we reduce the cost of the code that needs to
  // take image spacing into account.
  // It will be initialized later, once for each dimension.
  std::vector< SpacingType > divisionTable( maxSize, 0 );

  // Loop over all dimensions and compute the generalized distance transform
  // and voronoi map for each scanline.
  //
  // Scanning lines along dimension 0 exhibit a good hit rate of L1 cache and no
  // special care has to be taken to make use of this fact.
  //
  // We can work on multiple scanlines in parallel, however, thus making use of
  // multiple processors.
  // dimension 0
  const SpacingType sqrs = spacing[0]*spacing[0];
  Parabolas envelope;
  envelope.reserve(size[0]);

  // dimension 0, nice and easy
  //
  // To use openmp, we let distanceIt advance in direction 1 and parallelize
  // the loop over each scanline along direction 1. See below for details.
  distanceIt.SetDirection(1);
  distanceIt.GoToBegin();

  if (this->m_CreateVoronoiMap)
    {
    voronoiMapIt.SetDirection(1);
    voronoiMapIt.GoToBegin();
    }

  if (this->m_UseImageSpacing)
    {
    updateDivisionTable(divisionTable, spacing[0], size[0]);
    }

  // Distance in pixels from one scanline to the next
  const size_t stride = size[ 0 ];

  // Iterate over all scanlines. The outer while advances slices, the for-loop
  // below iterates along direction 1 over each scanline in direction 0.
  while (!distanceIt.IsAtEnd())
    {
    // Compute the generalized distance transform for the current scanline.
    // We are using raw pointers to the data in order to simplify access with
    // OpenMP.
    //
    // \todo When ITK offers new methods to store images in memory, this code
    // would need to be revisited, as it assumes row mayor layout and
    // contiguous memory.
    DistancePixelType *rawDistance = &distanceIt.Value();
    LabelPixelType *rawLabel =
      this->m_CreateVoronoiMap ?  &voronoiMapIt.Value() : 0;
    for (size_t i = 0; i < size[ 1 ]; ++i)
      {
      // Where does the scanline start in rawDistance?
      size_t offset = i * stride;

      // First compute the lower envelope of parabolas
      envelope.clear();
      if (this->m_UseImageSpacing)
        {
        for (size_t j = 0; j < size[0]; ++j, ++offset)
          {
          addParabola(
            envelope, size[0],
            j,
            *(rawDistance + offset),
            this->m_CreateVoronoiMap ?
            *(rawLabel + offset) :
            typename LabelImageType::PixelType(),
            divisionTable);
          }
        offset -= size[0];
        }
      else
        {
        for (size_t j = 0; j < size[ 0 ]; ++j, ++offset)
          {
          addParabola(
            envelope, size[0],
            j,
            *(rawDistance + offset),
            this->m_CreateVoronoiMap ?
            *(rawLabel + offset) :
            typename LabelImageType::PixelType());
          }
        offset -= size[0];
        }

      // And now sample the lower envelope for the whole scanline
      if (this->m_UseImageSpacing)
        {
        sampleValues(envelope, 0, size[0], rawDistance + offset, sqrs);
        }
      else
        {
        sampleValues(envelope, 0, size[0], rawDistance + offset);
        }

      if (this->m_CreateVoronoiMap)
        {
        sampleVoronoi(envelope, 0, size[0], rawLabel + offset);
        }
      }

    // Advance the iterators to the next slice
    for (size_t i = 0; i < size[ 1 ]; ++i)
      {
      distanceIt.NextLine();
      if (this->m_CreateVoronoiMap)
        {
        voronoiMapIt.NextLine();
        }
      }
    }


  // Dimensions 1 and up are interesting when it comes to L1 cache usage.
  //
  // We will build several envelopes in parallel along direction d by reading a
  // number of scanlines along direction 0 at once.
  //
  // For writing the result, we first sample all new scanlines into a cache
  // efficient buffer and then write that buffer into the output image, again in
  // a manner that is keeping L1 cache busy.
  //
  // This extra step is necessary because the loop for sampling the output can
  // be better optimized when it can be handled in one go. Therefore we can not
  // sample multiple envelopes at once without loosing a lot of runtime
  // performance.

  // The number of scanlines is the number of pixels that fit in a L1 cache
  // line or a small multiple thereof. This way, we will use everything from L1
  // cache (when the memory is properly aligned).
  const DistSizeValueType pixelsInCacheLine =
    std::max( (DistSizeValueType)1,
      (DistSizeValueType)(m_CacheLineSize/
        sizeof(typename DistanceImageType::PixelType)));

  // We are using a buffer for each thread for intermediate output.
  std::vector< DistancePixelType > valuesBuffer;
  valuesBuffer.resize(pixelsInCacheLine*maxSize);
  std::vector< LabelPixelType > voronoiBuffer;
  if( this->m_UseImageSpacing )
    {
    voronoiBuffer.resize(pixelsInCacheLine*maxSize);
    }

  // We are looping over the remaining image dimensions starting with the
  // highest. This is due to the following observation:
  //
  // The algorithm detects scanlines with all background pixels and can
  // afterwards skip writing them back. On the other hand, when you have a
  // single non-background pixel in an axial slice, after scanning and
  // resampling along dimension 0 you will have a whole scanline of
  // non-background pixels. After that, none of the scanlines on this slice in
  // dimension 1 will stay empty.
  //
  // The CT data sets that we have investigated almost never had a completely
  // empty axial slice. They had a certain amount of empty coronary slices,
  // however.
  //
  // This iteration order saves about 8% runtime on our test system, so we opted
  // to investigate that direction first.
  for (int d = FunctionImageType::ImageDimension-1; d > 0; --d)
  // Curious how it works for you the other way round? Then try this instead:
  //for (int d = 1; d < FunctionImageType::ImageDimension; ++d)
    {
    // Set up spacing and division info
    const SpacingType sqrsD = spacing[d]*spacing[d];

    if (this->m_UseImageSpacing)
      {
      updateDivisionTable(divisionTable, spacing[d], size[d]);
      }

    distanceIt.SetDirection(d);
    distanceIt.GoToBegin();

    if (this->m_CreateVoronoiMap)
      {
      voronoiMapIt.SetDirection(d);
      voronoiMapIt.GoToBegin();
      }

    std::vector< Parabolas > envelopeD( pixelsInCacheLine );
    for (size_t i = 0; i < pixelsInCacheLine; ++i)
      {
      envelopeD[i].reserve( maxSize );
      }

    // We call the group of scanlines that is sampled in parallel a strip. How
    // many strips are there? This many:
    const int strips =
      static_cast< int >(std::ceil(double(size[0]) / pixelsInCacheLine));

    // Compute stride from one scanline in direction d to the next
    size_t strideD = 1;
    for (int i = 0; i < d; ++i)
      {
      strideD *= size[i];
      }

    while (!distanceIt.IsAtEnd())
      {
      // Inside the for loop, we are working with lean and mean pointers into
      // the image buffer and offsets therein.
      // This makes parallelization with openmp easier as we do not need to
      // modify an iterator in the inner loop.
      //
      // \todo Again, this needs to be reinvestigated when other memory layouts
      // are introduced into ITK
      DistancePixelType *rawDistance = &distanceIt.Value();
      LabelPixelType *rawLabel =
        this->m_CreateVoronoiMap ?  &voronoiMapIt.Value() : 0;

      // When you compile with openmp support, this loop can be handled in
      // parallel. Each thread gets its own copy of envelope and output buffers
      // to work on. All other data can be shared as it is not written to.
      //
      // The alternative implementation to make valuesBuffer and voronoiBuffer
      // cover the whole slice and share the data for output has been discarded
      // as it would create a prohibitively large intermediate image when used
      // with large 2D-images.
      //
      // The static schedule seemed best for the task at hand. It has the least
      // overhead to implement and as each strip takes approximately the same
      // time to complete, it fits well.
      //
      // Each thread will handle a single strip and therefore have everything in
      // L1 cache of its CPU.
      for (int i = 0; i < strips; ++i)
        {
        const IndexValueType currentStripOffset = i * pixelsInCacheLine;

        // We can work on at most linesInCache lines in parallel to be cache
        // effective. To avoid wrap-around effects at the end of a line, we also
        // need to take the image size in dimension 0 in account.
        const DistSizeValueType parallelLines = std::min(
          pixelsInCacheLine,
          size[0] - i * pixelsInCacheLine);

        // Now compute the lower envelope of parabolas for each scanline
        for (size_t line = 0; line < parallelLines; ++line)
          {
          envelopeD[line].clear();
          }

        // Read a strip of parallelLines scanlines along direction d
        for (size_t j = 0, offset = currentStripOffset;
              j < size[d];
              ++j, offset += strideD - parallelLines)
          {
          if (this->m_UseImageSpacing)
            {
            // This is the inner loop where it matters to utilize L1 cache
            for (size_t line = 0; line < parallelLines; ++line, ++offset)
              {
              addParabola(
                envelopeD[line], size[d],
                j,
                *(rawDistance + offset),
                this->m_CreateVoronoiMap ?
                *(rawLabel + offset) :
                typename LabelImageType::PixelType(),
                divisionTable);
              }
            }
          else
            {
            // No wait, this one! ;-)
            for (size_t line = 0; line < parallelLines; ++line, ++offset)
              {
              addParabola(
                envelopeD[line], size[d],
                j,
                *(rawDistance + offset),
                this->m_CreateVoronoiMap ?
                *(rawLabel + offset) :
                typename LabelImageType::PixelType());
              }
            }
          }

        // And now evaluate the lower envelope
        //
        // We first sample into the intermediate buffer. Each envelope samples a
        // whole line.  This allows for efficient loop unrolling in sampleValues
        // by the compiler and effective use of L1 cache during the sampling.
        //
        // We can avoid to write empty lines into the output buffer because we
        // already know that the output image has to be all background anyway.
        // We only have to remember that a line is empty.
        //
        // There is also some possibility that we can avoid to write the whole
        // strip.
        //
        // In the output buffers, the scanlines are maxSize elements apart.
        bool stripNeedsCopy = true;
        std::vector< bool > lineNeedsCopy( parallelLines, false );
        if (m_UseImageSpacing)
          {
          for (size_t line = 0, offset = 0;
                line < parallelLines;
                ++line, offset += maxSize)
            {
            lineNeedsCopy[line] =
              sampleValues(
                envelopeD[line], 0, size[d], &valuesBuffer[offset], sqrsD);
            stripNeedsCopy &= lineNeedsCopy[line];
            }
          }
        else
          {
          for( size_t line = 0, offset = 0;
            line < parallelLines; ++line, offset += maxSize )
            {
            lineNeedsCopy[line] =
              sampleValues(
                envelopeD[line], 0, size[d], &valuesBuffer[offset]);
            stripNeedsCopy &= lineNeedsCopy[line];
            }
          }

        // Now we write the buffer to the output image. This is done in parallel
        // for all scanlines in the strip. Therefore, we utilize the cache lines
        // that hold a few words of each scanline in valuesBuffer and a cache
        // line for the output to rawDistance.
        //
        // With probably 1000 lines of L1 cache available and much less lines in
        // a single strip, we utilize L1 cache quite well again.
        if (!stripNeedsCopy) continue;
        for (size_t j = 0, outoffset = currentStripOffset;
          j < size[d]; ++j, outoffset += strideD - parallelLines)
          {
          for (size_t line = 0, inoffset = j;
            line < parallelLines;
            ++line,
            inoffset += maxSize, // ...scanlines are maxSize elements apart
            ++outoffset)
            {
            // outoffset advances by 1, so rawDistance is written sequentially.
            // inoffset advances by a larger amount, but after parallelLines
            // iterations, the cache line read first can be reused again.
            if (lineNeedsCopy[line])
              {
              *(rawDistance + outoffset) = valuesBuffer[inoffset];
              }
            }
          }

        if (this->m_CreateVoronoiMap)
          {
          // More of the same
          for (size_t line = 0, offset = 0;
            line < parallelLines; ++line, offset += maxSize)
               {
               if( lineNeedsCopy[line] )
                 {
                 sampleVoronoi(envelopeD[line], 0, size[d],
                   &voronoiBuffer[offset]);
                 }
               }

          for( size_t j = 0, outoffset = currentStripOffset;
            j < size[d];
            ++j, outoffset += strideD - parallelLines)
            {
            for( size_t line = 0, inoffset = j;
              line < parallelLines;
              ++line, inoffset += maxSize, ++outoffset)
              {
              if( lineNeedsCopy[line] )
                {
                *(rawLabel + outoffset) = voronoiBuffer[inoffset];
                }
              }
            }
          }
        }

      // Ok, all strips are done
      for (size_t i = 0; i < size[0]; ++i)
        {
        distanceIt.NextLine();
        if( this->m_CreateVoronoiMap )
          {
          voronoiMapIt.NextLine();
          }
        }
      }
    }
} // end GenerateData()

template < class TFunctionImage, class TDistanceImage, class TLabelImage >
typename GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::AbscissaIndexType
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::intersection(const Parabola &p, const Parabola &q,
  const std::vector< SpacingType > & divisionTable,
  const AbscissaIndexType &from,
  const AbscissaIndexType &to)
{
  // Parabolas must be different, because the intersection abscissa is not
  // defined otherwise
  assert(p.i != q.i);

  // We have fp(x) = (x - px)^2 + py and fq(x) = (x - qx)^2 + qy
  //     fp(x) = fq(x)
  // <=> (x - px)^2 + py = (x - qx)^2 + qy
  // <=> x^2 - 2xpx + px^2 + py = x^2 - 2xqx + qx^2 + qy
  // <=> -2xpx + 2xqx = qx^2 - px^2 + qy - py
  // <=> x(2qx - 2px) = qx^2 - px^2 + qy - py
  // <=> x = (qx^2 - px^2 + qy - py) / (2 (qx - px))
  // <=> x = 1/2 * ((qx^2 - px^2) / (qx - px) + (qy - py) / (qx - px))
  // <=> x = 1/2 * (qx + px + (qy - py) / (qx - px))
  //
  //   using the index i and spacing s instead of x (x = i s)
  //
  // <=> i = 1/2 * (qi + pi + (qy - py) / (s^2(qi - pi)))
  //
  // This has to be evaluated in floating point precision, because
  // the denominator of the fraction can be < 1.0 for small spacings.
  //
  // We return the largest index where p is strictly below q, therefore we
  // truncate and increment
  //  i = floor(1/2 * (qi + pi + (qy - py) / (s^2(qi - pi)))) + 1
  //
  // Furthermore, i will be clamped to the range [from, to].

  // Shortcut for same y:
  //
  // Return the mean of the two indices.
  if (p.y == q.y)
    {
    return (p.i + q.i) / 2 + 1;
    }

  // Estimation whether i <= from.
  //
  // <=> floor(1/2 * (qi + pi + (qy - py) / (s^2(qi - pi)))) + 1 <= from
  // <=> floor(1/2 * (qi + pi + (qy - py) / (s^2(qi - pi)))) <= from - 1
  // <=  1/2 * (qi + pi + (qy - py) / (s^2(qi - pi))) <= from - 1
  // <=> qi + pi + (qy - py) / (s^2(qi - pi)) <= 2*from - 2
  // <=  2*qi + (qy - py) / (s^2(qi - pi)) <= 2*from - 2
  //     and for non-positive qy-py:
  // <=  2*qi <= 2*from - 2
  // <=> qi <= from - 1
  // <=> qi <  from
  //
  // However, the optimization is a bad idea, because the extra logic costs
  // more runtime than it saves.
  //  if (qy-py <= 0 && q.i < from)
  //    return from;
  const AbscissaIndexType i(
      (static_cast<SpacingType>(q.i+p.i) +
       divisionTable[ q.i-p.i ] * static_cast<SpacingType>(q.y-p.y)) *
      0.5 + 1);
  if( i < from )
    {
    return from;
    }
  else if( i > to )
    {
    return to;
    }
  else
    {
    return i;
    }
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
typename GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::AbscissaIndexType
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::intersection(const Parabola &p, const Parabola &q,
             const AbscissaIndexType &from,
             const AbscissaIndexType &to)
{
  // All comments from the function above apply, but this one does not use
  // spacing.
  assert(p.i != q.i);
  if (p.y == q.y)
    {
    return (p.i + q.i) / 2 + 1;
    }

  const AbscissaIndexType di = q.i - p.i;
  if (di != 1)
    {
    const AbscissaIndexType i = ((q.i + p.i) + (q.y - p.y) / di) / 2 + 1;
    if( i < from )
      {
      return from;
      }
    else if( i > to )
      {
      return to;
      }
    else
      {
      return i;
      }
    }
  else
    {
    const AbscissaIndexType i = ((q.i + p.i) + (q.y - p.y)) / 2 + 1;
    if( i < from )
      {
      return from;
      }
    else if( i > to )
      {
      return to;
      }
    else
      {
      return i;
      }
    }
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
void
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::addParabola(
  Parabolas &envelope,
  const AbscissaIndexType& to,
  const AbscissaIndexType& pi, const DistancePixelType& py,
  const LabelPixelType& pl,
  const std::vector< SpacingType > & divisionTable )
{
  // We do not care for background values
  if (py >= m_MaximalSquaredDistance)
    {
    return;
    }

  Parabola p = {pi, py, pl, 0};

  // Parabolas have to be added with increasing abscissas
  assert(envelope.empty() || p.i > envelope.back().i);

  // The region where the newly added parabola is minimal is extending to
  // infinity, because with the highest abscissa it will eventually be below any
  // of the other parabolas in the envelope.  We are finding the abscissa where
  // the minimal region is begining by intersection with the parabolas in the
  // envelope.
  AbscissaIndexType i = 0;
  while (!envelope.empty() &&
    (i = intersection( envelope.back(), p, divisionTable,
      envelope.back().dominantFrom, to)) <= envelope.back().dominantFrom)
    {
    // The new parabola is below the whole last parabola region. Hence, the
    // last parabola region is not part of the envelope anymore.
    envelope.pop_back();
    }

  // All parabolas that are made obsolete by the new one have been removed
  // and the new one can be inserted.
  //
  // We are only interested in intersections to the left of 'to'.
  // Parabolas to the right of it will not be sampled anyway.
  if (i < to)
    {
    p.dominantFrom = i;
    envelope.push_back(p);
    }
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
void
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::addParabola(
  Parabolas &envelope,
  const AbscissaIndexType& to,
  const AbscissaIndexType& pi, const DistancePixelType& py,
  const LabelPixelType& pl)
{
  // More of the same, without spacing.
  if (py >= m_MaximalSquaredDistance)
    {
    return;
    }

  Parabola p = {pi, py, pl, 0};
  assert(envelope.empty() || p.i > envelope.back().i);
  AbscissaIndexType i = 0;
  while (!envelope.empty() &&
    (i = intersection( envelope.back(), p,
      envelope.back().dominantFrom, to)) <= envelope.back().dominantFrom)
    {
    envelope.pop_back();
    }

  if (i < to)
    {
    p.dominantFrom = i;
    envelope.push_back(p);
    }
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
bool
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::sampleValues(
  const Parabolas &envelope,
  const AbscissaIndexType &from, const long &steps,
  DistancePixelType *buffer,
  const SpacingType &sqrs)
{
  // An empty envelope signifies that only background values are present in
  // the scanline. Nothing needs to be sampled at all
  if (envelope.empty())
    {
    return false;
    }

  // Search for the correct parabola
  typename Parabolas::size_type parabolaIndex = 0;
  while (parabolaIndex != envelope.size() - 1 &&
    envelope[parabolaIndex+1].dominantFrom <= from)
    {
    ++parabolaIndex;
    }

  // To avoid casting in the following inner loop, we set up a variable 'delta'
  // with SpacingType that should hold integers in the range [0,steps).
  // Therefore, we rely on the fact that it can hold all expected values
  // accurately.
  //
  // This is the case for float/double and scanlines with a size of up to
  // 2^23 pixels, but just to be sure:
  assert(pow(SpacingType(10),
    std::numeric_limits<SpacingType>::digits10) >= steps);
  // If the assert fails, you need to make delta an instance of
  // AbscissaIndexType at the cost of an added cast in the inner loop.

  // Loop over all parabolas covering the interval [from, from+steps)
  for (AbscissaIndexType i = from;
    i != from + steps; ++parabolaIndex) // i is increased within the loop
    {
    // Compute the interval for which the current parabola needs to be sampled.
    // This goes from i either to the start of the dominance interval of the
    // next parabola or to the end of the scanline.
    const Parabola &p = envelope[parabolaIndex];

    AbscissaIndexType stepsInRegion;
    if( parabolaIndex < envelope.size() - 1 &&
      envelope[parabolaIndex+1].dominantFrom <= from+steps )
      {
      stepsInRegion = envelope[parabolaIndex+1].dominantFrom - i;
      }
    else
      {
      stepsInRegion = from+steps - i;
      }

    // addParabola needs to make sure that every parabola has at least one
    // pixel where it is relevant.
    assert(stepsInRegion > 0);

    // Initialize sampling
    SpacingType delta = i-p.i;

    // A consecutive number of foreground pixels lead to a consecutive
    // number of parabolas, each with a dominance region that only covers
    // their apex, i.e. p.i. In that case, we do not need to sample the
    // parabola as we already know its value at this point.
    //
    // This happens often enough to justify the shortcut.
    if (delta == 0 && stepsInRegion == 1)
      {
      *buffer++ = p.y;
      ++i;
      continue;
      }

    // Otherwise we sample correctly
    const SpacingType py(p.y);
    i += stepsInRegion; // We don't need i any more, increase here.
    for (; stepsInRegion; --stepsInRegion, ++buffer, ++delta)
      {
      DistancePixelType value(py + delta * delta * sqrs);
      if( value < m_MaximalSquaredDistance )
        {
        *buffer = value;
        }
      else
        {
        *buffer = m_MaximalSquaredDistance;
        }
      }
    }
  return true;
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
bool
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::sampleValues(
  const Parabolas &envelope,
  const AbscissaIndexType &from, const long &steps,
  DistancePixelType *buffer)
{
  // More of the same, but see below for sampling in the integer domain
  if (envelope.empty())
    {
    return false;
    }

  typename Parabolas::size_type parabolaIndex = 0;

  while (parabolaIndex != envelope.size() - 1 &&
    envelope[parabolaIndex+1].dominantFrom <= from)
    {
    ++parabolaIndex;
    }

  for (AbscissaIndexType i = from;
    i != from + steps; ++parabolaIndex)
    {
    const Parabola &p = envelope[parabolaIndex];
    AbscissaIndexType stepsInRegion =
      (parabolaIndex < envelope.size() - 1 &&
       envelope[parabolaIndex+1].dominantFrom <= from+steps) ?
      envelope[parabolaIndex+1].dominantFrom - i :
      from+steps - i;

    AbscissaIndexType delta = i-p.i;
    if (delta == 0 && stepsInRegion == 1)
      {
      *buffer++ = p.y;
      ++i;
      continue;
      }

    i += stepsInRegion;
    DistancePixelType value = delta*delta + p.y; // Integer only, we can add p.y right now
    AbscissaIndexType inc = 2*delta+1;
    for (; stepsInRegion; --stepsInRegion, ++buffer)
      {
      // p.y is preadded, spacing is 1
      if( value < m_MaximalSquaredDistance )
        {
        *buffer = value;
        }
      else
        {
        *buffer = m_MaximalSquaredDistance;
        }

      // Advance the sample position, update value to current height
      // (x+1)^2 + b = (x^2+2x+1) + b.
      // We started with x == delta and set up inc in order to do the
      // sampling with two additions instead of one multiplication
      value += inc;
      inc += 2;
      }
    }
  return true;
}

template < class TFunctionImage,class TDistanceImage, class TLabelImage >
bool
GeneralizedDistanceTransformImageFilter<
  TFunctionImage, TDistanceImage, TLabelImage >
::sampleVoronoi(
  const Parabolas &envelope,
  const AbscissaIndexType &from, const long &steps,
  LabelPixelType *buffer)
{
  if (envelope.empty())
    {
    return false;
    }

  typename Parabolas::size_type parabolaIndex = 0;
  // Search for the correct parabola
  // We are only interested in the labels, therefore we can join parabola
  // dominance intervals that represent the same region.
  while (parabolaIndex != envelope.size() - 1 &&
    (envelope[parabolaIndex+1].dominantFrom <= from ||
    envelope[parabolaIndex+1].l == envelope[parabolaIndex].l))
    {
    ++parabolaIndex;
    }

  AbscissaIndexType i = from;
  while (i != from + steps)
    {
    const Parabola &p = envelope[parabolaIndex];

    AbscissaIndexType stepsInRegion;

    if( parabolaIndex < envelope.size() - 1 &&
      envelope[parabolaIndex+1].dominantFrom <= from+steps )
      {
      stepsInRegion = envelope[parabolaIndex+1].dominantFrom - i;
      }
    else
      {
      stepsInRegion = from+steps - i;
      }

    std::fill_n(buffer, stepsInRegion, p.l);
    i += stepsInRegion;
    buffer += stepsInRegion;

    // Increase the parabola index to the next parabola of a different label
    ++parabolaIndex;
    while (parabolaIndex < envelope.size() - 1 &&
           envelope[parabolaIndex+1].l == envelope[parabolaIndex].l)
      {
      ++parabolaIndex;
      }
    }
  return true;
}
} // end namespace itk
#endif
