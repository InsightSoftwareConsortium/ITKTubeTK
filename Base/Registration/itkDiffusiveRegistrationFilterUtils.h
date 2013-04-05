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
#ifndef __itkDiffusiveRegistrationFilterUtils_h
#define __itkDiffusiveRegistrationFilterUtils_h

namespace itk
{

class ITK_EXPORT DiffusiveRegistrationFilterUtils
{
public:

  /** Helper function to allocate an image based on a template */
  template< class UnallocatedImagePointer, class TemplateImagePointer >
  static void AllocateSpaceForImage( UnallocatedImagePointer & image,
                                     const TemplateImagePointer & templateImage );

  /** Helper function to check whether the attributes of an image match a
    * template */
  template< class CheckedImageType, class TemplateImageType >
  static bool CompareImageAttributes( const CheckedImageType * image,
                                      const TemplateImageType * templateImage );

  /** Resamples an image to a template using nearest neighbor interpolation */
  template< class ResampleImagePointer, class TemplateImagePointer >
  static void ResampleImageNearestNeighbor(
      const ResampleImagePointer & highResolutionImage,
      const TemplateImagePointer & templateImage,
      ResampleImagePointer & resampledImage );

  /** Resamples an image to a template using linear interpolation */
  template< class ResampleImagePointer, class TemplateImagePointer >
  static void ResampleImageLinear(
      const ResampleImagePointer & highResolutionImage,
      const TemplateImagePointer & templateImage,
      ResampleImagePointer & resampledImage );

  /** Resamples a vector image to a template using linear interpolation.  If
   *  normalize is true, the vectors will be scaled to length 1 after the
   *  resampling. */
  template< class VectorResampleImagePointer, class TemplateImagePointer >
  static void VectorResampleImageLinear(
      const VectorResampleImagePointer & highResolutionImage,
      const TemplateImagePointer & templateImage,
      VectorResampleImagePointer & resampledImage,
      bool normalize = false );

  /** Normalizes a vector field to ensure each vector has length 1 */
  template< class VectorImagePointer >
  static void NormalizeVectorField( VectorImagePointer & image );

  /** Computes the minimum and maximum intensity in an image */
  template< class ImageType >
  static bool IsIntensityRangeBetween0And1( ImageType * image );

  /** Extracts the x, y, z components of a deformation field. */
  template< class DeformationFieldType, class DeformationComponentImageArrayType >
  static void ExtractXYZComponentsFromDeformationField(
      const DeformationFieldType * deformationField,
      DeformationComponentImageArrayType & deformationComponentImages );

};


/** Struct to simply get the face list and an iterator over the face list
 *  when processing an image.  Designed for use with SmartPointers. */
template< class ImageType >
struct FaceStruct
{
  typedef NeighborhoodAlgorithm::ImageBoundaryFacesCalculator
      < typename ImageType::ObjectType >                  FaceCalculatorType;
        typedef typename FaceCalculatorType::FaceListType FaceListType;
        typedef typename FaceListType::iterator           FaceListIteratorType;

  FaceStruct()
    {
    numberOfTerms = 0;
    }

  FaceStruct( const ImageType& image,
              typename ImageType::ObjectType::RegionType region,
              typename ImageType::ObjectType::SizeType radius )
    {
    numberOfTerms = 0;
    if( image.GetPointer() )
      {
      faceLists.push_back( faceCalculator( image, region, radius ) );
      numberOfTerms = 1;
      }
    }

  FaceStruct( const std::vector< ImageType >& images,
              typename ImageType::ObjectType::RegionType region,
              typename ImageType::ObjectType::SizeType radius )
    {
    numberOfTerms = 0;
    for( int i = 0; i < (int) images.size(); i++ )
      {
      if( images[i].GetPointer() )
        {
        faceLists.push_back( faceCalculator( images[i], region, radius ) );
        numberOfTerms++;
        }
      }
    }

  template< unsigned int VLength >
  FaceStruct( const itk::FixedArray< ImageType, VLength >& images,
              typename ImageType::ObjectType::RegionType region,
              typename ImageType::ObjectType::SizeType radius )
    {
    numberOfTerms = 0;
    for( int i = 0; i < (int) images.Size(); i++ )
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
              std::vector< itk::FixedArray< ImageType, VLength > > &images,
              typename ImageType::ObjectType::RegionType region,
              typename ImageType::ObjectType::SizeType radius )
    {
    numberOfTerms = 0;
    for( int i = 0; i < (int) images.size(); i++)
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

  void GoToBegin()
    {
    if( (int) faceListIts.size() != numberOfTerms )
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

  bool IsAtEnd()
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

  void Increment()
    {
    for( int i = 0; i < numberOfTerms; i++ )
      {
      ++faceListIts[i];
      }
    }

  template< class IteratorType >
  void SetIteratorToCurrentFace(
      IteratorType& iterator,
      const ImageType& image,
      typename ImageType::ObjectType::SizeType radius )
    {
    if( image.GetPointer() )
      {
      iterator = IteratorType( radius, image, *faceListIts[0] );
      }
    else
      {
      iterator = IteratorType();
      }
    }

  template< class IteratorType >
  void SetIteratorToCurrentFace(
      IteratorType& iterator,
      const ImageType& image )
    {
    if( image.GetPointer() )
      {
      iterator = IteratorType( image, *faceListIts[0] );
      }
    else
      {
      iterator = IteratorType();
      }
    }

  template< class IteratorType >
  void SetIteratorToCurrentFace(
      std::vector< IteratorType >& iterators,
      const std::vector< ImageType >& images,
      typename ImageType::ObjectType::SizeType radius )
    {
    if( (int) iterators.size() != numberOfTerms )
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        if( images[i].GetPointer() )
          {
          iterators.push_back(
              IteratorType( radius, images[i], *faceListIts[i] ) );
          }
        else
          {
          iterators.push_back( IteratorType() );
          }
        }
      }
    else
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        if( images[i].GetPointer() )
          {
          iterators[i] = IteratorType( radius, images[i], *faceListIts[i] );
          }
        else
          {
          iterators[i] = IteratorType();
          }
        }
      }
    }

  template< class IteratorType >
  void SetIteratorToCurrentFace(
      std::vector< IteratorType >& iterators,
      const std::vector< ImageType >& images )
    {
    if( (int) iterators.size() != numberOfTerms )
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        if( images[i].GetPointer() )
          {
          iterators.push_back( IteratorType( images[i], *faceListIts[i] ) );
          }
        else
          {
          iterators.push_back( IteratorType() );
          }
        }
      }
    else
      {
      for( int i = 0; i < numberOfTerms; i++ )
        {
        if( images[i].GetPointer() )
          {
          iterators[i] = IteratorType( images[i], *faceListIts[i] );
          }
        else
          {
          iterators[i] = IteratorType();
          }
        }
      }
    }

    template< class IteratorType, unsigned int VLength >
    void SetIteratorToCurrentFace(
        std::vector< itk::FixedArray< IteratorType, VLength > > &iterators,
        const std::vector< itk::FixedArray< ImageType, VLength > > & images )
    {
    int c = 0;
    if( (int) iterators.size() != (int) images.size() )
        {
        for( int i = 0; i < (int) images.size(); i++ )
          {
          itk::FixedArray< IteratorType, VLength > fixedArray;
          for( int j = 0; j < (int) images[i].Size(); j++ )
            {
            if( images[i][j] )
              {
              fixedArray[j] = IteratorType( images[i][j], *faceListIts[c] );
              c++;
              }
            else
              {
              fixedArray[j] = IteratorType();
              }
            }
          iterators.push_back( fixedArray );
          }
        }
      else
        {
        for( int i = 0; i < (int) images.size(); i++ )
          {
          itk::FixedArray< IteratorType, VLength > fixedArray;
          for( int j = 0; j < (int) images[i].Size(); j++ )
            {
            if( images[i][j].GetPointer() )
              {
              fixedArray[j] = IteratorType( images[i][j], *faceListIts[c] );
              c++;
              }
            else
              {
              fixedArray[j] = IteratorType();
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
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
# include "itkDiffusiveRegistrationFilterUtils.txx"
#endif

#endif
