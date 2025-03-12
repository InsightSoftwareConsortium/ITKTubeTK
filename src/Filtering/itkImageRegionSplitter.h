/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkImageRegionSplitter_h
#define itkImageRegionSplitter_h

#include "itkImageRegion.h"
#include "itkImageRegionSplitterBase.h"

namespace itk
{
/** \class ImageRegionSplitter
 * \brief Divide a region into several pieces.
 *
 * ImageRegionSplitter divides an ImageRegion into smaller regions.
 * ImageRegionSplitter is used by the StreamingImageFilter to divide a
 * requested output region into a series of smaller requests of the
 * pipeline.  This object has two basic methods: GetNumberOfSplits()
 * and GetSplit().
 *
 * GetNumberOfSplits() is used to determine how may subregions a given
 * region can be divided.  You call GetNumberOfSplits with an argument
 * that is the number of subregions you want.  If the image region can
 * support that number of subregions, that number is returned.
 * Otherwise, the maximum number of splits less then or equal to the
 * argumen  be returned.  For example, if a region splitter class only divides
 * a region into horizontal slabs, then the maximum number of splits
 * will be the number of rows in the region.
 *
 * GetSplit() returns the ith of N subregions ( as an ImageRegion object ).
 *
 * This ImageRegionSplitter class divides a region along the outermost
 * dimension. If the outermost dimension has size 1 ( i.e. a volume
 * with a single slice ), the ImageRegionSplitter will divide the
 * region along the next outermost dimension. If that dimension has size 1,
 * the process continues with the next outermost dimension.
 *
 * Other ImageRegionSplitter subclasses could divide an image into
 * more uniform shaped regions instead of slabs.
 *
 * \deprecated The new class ImageRegionSplitterSlowDimension can be
 * used as a drop in replacement for functionality as it
 * implements the same algorithm. The ImageRegionSplitterBase is now
 * the abstract base class for all image splitter. New splitting
 * object should be derived from that class.
 *
 * \sa ImageRegionSplitterSlowDimension
 * \sa ImageRegionSplitterBase
 *
 * \ingroup ITKDeprecated
 */

template <unsigned int VImageDimension>
class ImageRegionSplitter : public ImageRegionSplitterBase
{
public:
  /** Standard class type alias. */
  using Self = ImageRegionSplitter;
  using Superclass = ImageRegionSplitterBase;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro(ImageRegionSplitter);

  /** Dimension of the image available at compile time. */
  static constexpr unsigned int ImageDimension = VImageDimension;

  /** Dimension of the image available at run time. */
  static unsigned int
  GetImageDimension()
  {
    return VImageDimension;
  }

  /** Index type alias support. An index is used to access pixel values. */
  using IndexType = Index<VImageDimension>;

  /** Size type alias support. A size is used to define region bounds. */
  using SizeType = Size<VImageDimension>;
  using SizeValueType = typename SizeType::SizeValueType;

  /** Region type alias support.   */
  using RegionType = ImageRegion<VImageDimension>;

  /** How many pieces can the specifed region be split? A given region
   * cannot always be divided into the requested number of pieces.  For
   * instance, if the numberOfPieces exceeds the number of pixels along
   * a certain dimensions, then some splits will not be possible. This
   * method returns a number less than or equal to the requested number
   * of pieces. */
  virtual unsigned int
  GetNumberOfSplits(const RegionType & region, unsigned int requestedNumber);


  /** Get a region definition that represents the ith piece a specified region.
   * The "numberOfPieces" must be equal to what
   * GetNumberOfSplits() returns. */
  virtual RegionType
  GetSplit(unsigned int i, unsigned int numberOfPieces, const RegionType & region);

protected:
  ImageRegionSplitter() {}
  ~ImageRegionSplitter() {}

  virtual unsigned int
  GetNumberOfSplitsInternal(unsigned int,
                            const IndexValueType regionIndex[],
                            const SizeValueType  regionSize[],
                            unsigned int         requestedNumber) const override
  {
    // this function adapts the legecy method, defined in this class
    // be used by the ImageRegionSplitterBase.
    IndexType idx;
    idx.SetIndex(regionIndex);
    SizeType sz;
    sz.SetSize(regionSize);
    RegionType region = RegionType(idx, sz);

    Self * nonconst_this = const_cast<Self *>(this);
    return nonconst_this->GetNumberOfSplits(region, requestedNumber);
  }

  virtual unsigned int
  GetSplitInternal(unsigned int   dim,
                   unsigned int   i,
                   unsigned int   numberOfPieces,
                   IndexValueType regionIndex[],
                   SizeValueType  regionSize[]) const override
  {
    // this function adapts the legecy method, defined in this class
    // be used by the ImageRegionSplitterBase.
    IndexType idx;
    idx.SetIndex(regionIndex);
    SizeType sz;
    sz.SetSize(regionSize);
    RegionType region = RegionType(idx, sz);

    Self * nonconst_this = const_cast<Self *>(this);
    region = nonconst_this->GetSplit(i, numberOfPieces, region);

    for (unsigned int d = 0; d < dim; ++d)
    {
      regionIndex[d] = region.GetIndex(d);
      regionSize[d] = region.GetSize(d);
    }
    return numberOfPieces;
  }

  virtual void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  // purposely not implemented
  ImageRegionSplitter(const ImageRegionSplitter &);

  // purposely not implemented
  void
  operator=(const ImageRegionSplitter &);
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkImageRegionSplitter.hxx"
#endif

#endif
