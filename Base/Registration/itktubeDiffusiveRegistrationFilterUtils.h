/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeDiffusiveRegistrationFilterUtils_h
#define __itktubeDiffusiveRegistrationFilterUtils_h

#include <vector>

namespace itk
{

namespace tube
{

class DiffusiveRegistrationFilterUtils
{
public:

  /** Helper function to allocate an image based on a template */
  template< class TUnallocatedImagePointer, class TTemplateImagePointer >
  static void AllocateSpaceForImage( TUnallocatedImagePointer & image,
                                     const TTemplateImagePointer & templateImage );

  /** Helper function to check whether the attributes of an image match a
    * template */
  template< class TCheckedImage, class TTemplateImage >
  static bool CompareImageAttributes( const TCheckedImage * image,
                                      const TTemplateImage * templateImage );

  /** Resamples an image to a template using nearest neighbor interpolation */
  template< class TResampleImagePointer, class TTemplateImagePointer >
  static void ResampleImageNearestNeighbor(
      const TResampleImagePointer & highResolutionImage,
      const TTemplateImagePointer & templateImage,
      TResampleImagePointer & resampledImage );

  /** Resamples an image to a template using linear interpolation */
  template< class TResampleImagePointer, class TTemplateImagePointer >
  static void ResampleImageLinear(
      const TResampleImagePointer & highResolutionImage,
      const TTemplateImagePointer & templateImage,
      TResampleImagePointer & resampledImage );

  /** Resamples a vector image to a template using linear interpolation.  If
   *  normalize is true, the vectors will be scaled to length 1 after the
   *  resampling. */
  template< class TVectorResampleImagePointer, class TTemplateImagePointer >
  static void VectorResampleImageLinear(
      const TVectorResampleImagePointer & highResolutionImage,
      const TTemplateImagePointer & templateImage,
      TVectorResampleImagePointer & resampledImage,
      bool normalize = false );

  /** Normalizes a vector field to ensure each vector has length 1 */
  template< class TVectorImagePointer >
  static void NormalizeVectorField( TVectorImagePointer & image );

  /** Computes the minimum and maximum intensity in an image */
  template< class TImage >
  static bool IsIntensityRangeBetween0And1( TImage * image );

  /** Extracts the x, y, z components of a deformation field. */
  template< class TDeformationField, class TDeformationComponentImageArray >
  static void ExtractXYZComponentsFromDeformationField(
      const TDeformationField * deformationField,
      TDeformationComponentImageArray & deformationComponentImages );

}; // End class DiffusiveRegistrationFilterUtils


/** Struct to simply get the face list and an iterator over the face list
 *  when processing an image.  Designed for use with SmartPointers. */
template< class TImage >
struct FaceStruct
{
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
      < typename TImage::ObjectType >                     FaceCalculatorType;
        typedef typename FaceCalculatorType::FaceListType FaceListType;
        typedef typename FaceListType::iterator           FaceListIteratorType;

  FaceStruct( void )
    {
    numberOfTerms = 0;
    }

  FaceStruct( const TImage& image,
              typename TImage::ObjectType::RegionType region,
              typename TImage::ObjectType::SizeType radius )
    {
    numberOfTerms = 0;
    if( image.GetPointer() )
      {
      faceLists.push_back( faceCalculator( image, region, radius ) );
      numberOfTerms = 1;
      }
    }

  FaceStruct( const std::vector< TImage >& images,
              typename TImage::ObjectType::RegionType region,
              typename TImage::ObjectType::SizeType radius )
    {
    numberOfTerms = 0;
    for( int i = 0; i < ( int ) images.size(); i++ )
      {
      if( images[i].GetPointer() )
        {
        faceLists.push_back( faceCalculator( images[i], region, radius ) );
        numberOfTerms++;
        }
      }
    }

  template< unsigned int VLength >
  FaceStruct( const itk::FixedArray< TImage, VLength >& images,
              typename TImage::ObjectType::RegionType region,
              typename TImage::ObjectType::SizeType radius )
    {
    numberOfTerms = 0;
    for( int i = 0; i < ( int ) images.Size(); i++ )
      {
      if( images[i].GetPointer() )
        {
        faceLists.push_back( faceCalculator( images[i], region, radius ) );
        numberOfTerms++;
        }
      }
    }

  template< unsigned int VLength >
  FaceStruct( const
              std::vector< itk::FixedArray< TImage, VLength > > &images,
              typename TImage::ObjectType::RegionType region,
              typename TImage::ObjectType::SizeType radius )
    {
    numberOfTerms = 0;
    for( int i = 0; i < ( int ) images.size(); i++ )
      {
      for( unsigned int j = 0; j < images[i].Size(); j++ )
        {
        if( images[i][j].GetPointer() )
          {
          faceLists.push_back( faceCalculator( images[i][j], region, radius ) );
          numberOfTerms++;
          }
        }
      }
    }

  void GoToBegin( void )
    {
    if( ( int ) faceListIts.size() != numberOfTerms )
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        faceListIts.push_back( faceLists[i].begin() );
        }
      }
    else
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        faceListIts[i] = faceLists[i].begin();
        }
      }
    }

  bool IsAtEnd( void )
    {
    for( int i = 0; i < numberOfTerms; i++ )
      {
      if( faceListIts[i] == faceLists[i].end() )
        {
        return true;
        }
      }
    return false;
    }

  void Increment( void )
    {
    for( int i = 0; i < numberOfTerms; i++ )
      {
      ++faceListIts[i];
      }
    }

  template< class TIterator >
  void SetIteratorToCurrentFace(
      TIterator& iterator,
      const TImage& image,
      typename TImage::ObjectType::SizeType radius )
    {
    if( image.GetPointer() )
      {
      iterator = TIterator( radius, image, *faceListIts[0] );
      }
    else
      {
      iterator = TIterator();
      }
    }

  template< class TIterator >
  void SetIteratorToCurrentFace(
      TIterator& iterator,
      const TImage& image )
    {
    if( image.GetPointer() )
      {
      iterator = TIterator( image, *faceListIts[0] );
      }
    else
      {
      iterator = TIterator();
      }
    }

  template< class TIterator >
  void SetIteratorToCurrentFace(
      std::vector< TIterator >& iterators,
      const std::vector< TImage >& images,
      typename TImage::ObjectType::SizeType radius )
    {
    if( ( int ) iterators.size() != numberOfTerms )
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        if( images[i].GetPointer() )
          {
          iterators.push_back(
              TIterator( radius, images[i], *faceListIts[i] ) );
          }
        else
          {
          iterators.push_back( TIterator() );
          }
        }
      }
    else
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        if( images[i].GetPointer() )
          {
          iterators[i] = TIterator( radius, images[i], *faceListIts[i] );
          }
        else
          {
          iterators[i] = TIterator();
          }
        }
      }
    }

  template< class TIterator >
  void SetIteratorToCurrentFace(
      std::vector< TIterator >& iterators,
      const std::vector< TImage >& images )
    {
    if( ( int ) iterators.size() != numberOfTerms )
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        if( images[i].GetPointer() )
          {
          iterators.push_back( TIterator( images[i], *faceListIts[i] ) );
          }
        else
          {
          iterators.push_back( TIterator() );
          }
        }
      }
    else
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        if( images[i].GetPointer() )
          {
          iterators[i] = TIterator( images[i], *faceListIts[i] );
          }
        else
          {
          iterators[i] = TIterator();
          }
        }
      }
    }

    template< class TIterator, unsigned int VLength >
    void SetIteratorToCurrentFace(
        std::vector< itk::FixedArray< TIterator, VLength > > &iterators,
        const std::vector< itk::FixedArray< TImage, VLength > > & images )
    {
    int c = 0;
    if( ( int ) iterators.size() != ( int ) images.size() )
        {
        for( int i = 0; i < ( int ) images.size(); i++ )
          {
          itk::FixedArray< TIterator, VLength > fixedArray;
          for( int j = 0; j < ( int ) images[i].Size(); j++ )
            {
            if( images[i][j] )
              {
              fixedArray[j] = TIterator( images[i][j], *faceListIts[c] );
              c++;
              }
            else
              {
              fixedArray[j] = TIterator();
              }
            }
          iterators.push_back( fixedArray );
          }
        }
      else
        {
        for( int i = 0; i < ( int ) images.size(); i++ )
          {
          itk::FixedArray< TIterator, VLength > fixedArray;
          for( int j = 0; j < ( int ) images[i].Size(); j++ )
            {
            if( images[i][j].GetPointer() )
              {
              fixedArray[j] = TIterator( images[i][j], *faceListIts[c] );
              c++;
              }
            else
              {
              fixedArray[j] = TIterator();
              }
            }
          iterators[i] = fixedArray;
          }
        }
      }

  FaceCalculatorType                     faceCalculator;
  std::vector< FaceListType >            faceLists;
  std::vector< FaceListIteratorType >    faceListIts;
  int                                    numberOfTerms;

}; // End struct FaceStruct

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeDiffusiveRegistrationFilterUtils.hxx"
#endif

#endif // End !defined( __itktubeDiffusiveRegistrationFilterUtils_h )
