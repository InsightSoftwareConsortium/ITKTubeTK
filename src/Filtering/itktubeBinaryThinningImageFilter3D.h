/*=========================================================================

Library:   TubeTKLib

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
#ifndef __itktubeBinaryThinningImageFilter3D_h
#define __itktubeBinaryThinningImageFilter3D_h

#include <list>
#include <itkNeighborhoodIterator.h>
#include <itkImageToImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkConstantBoundaryCondition.h>

namespace itk {

namespace tube {
/** \class BinaryThinningImageFilter3D
*
* \brief This filter computes one-pixel-wide skeleton of a 3D input image.
*
* This class is parametrized over the type of the input image
* and the type of the output image.
* 
* The input is assumed to be a binary image. All non-zero valued voxels
* are set to 1 internally to simplify the computation. The filter will
* produce a skeleton of the object.  The output background values are 0,
* and the foreground values are 1.
* 
* A 26-neighbourhood configuration is used for the foreground and a
* 6-neighbourhood configuration for the background. Thinning is performed
* symmetrically in order to guarantee that the skeleton lies medial within
* the object.
*
* This filter is a parallel thinning algorithm and is an implementation
* of the algorithm described in:
* 
* T.C. Lee, R.L. Kashyap, and C.N. Chu.
* Building skeleton models via 3-D medial surface/axis thinning algorithms.
* Computer Vision, Graphics, and Image Processing, 56(6):462--478, 1994.
* 
* To do: Make use of multi-threading.
*
* \author Hanno Homann, Oxford University, Wolfson Medical Vision Lab, UK.
* 
* Originally available via Insight Journal: https://hdl.handle.net/1926/1292
* and made available under the Creative Commons Attirubtion License v3.0
* https://creativecommons.org/licenses/by/3.0/legalcode
*
* Adapted for inclusion in TubeTK by Jared Vicory, Kitware, Inc., 18/07/2018
* 
* \sa MorphologyImageFilter
* \ingroup ImageEnhancement MathematicalMorphologyImageFilters
*/

template <class TInputImage, class TOutputImage=TInputImage>
class ITK_EXPORT BinaryThinningImageFilter3D :
    public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef BinaryThinningImageFilter3D                  Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  typedef typename TInputImage::PointType              PointType;

  typedef std::vector< PointType >                     EndPointListType;

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro( BinaryThinningImageFilter3D);

  /** Type for input image. */
  typedef   TInputImage       InputImageType;

  /** Type for output image: Skelenton of the object.  */
  typedef   TOutputImage      OutputImageType;

  /** Type for the region of the input image. */
  typedef typename InputImageType::RegionType RegionType;

  /** Type for the pixel type of the input image. */
  typedef typename InputImageType::PixelType InputImagePixelType ;

  /** Type for the pixel type of the input image. */
  typedef typename OutputImageType::PixelType OutputImagePixelType ;

  /** Type for the size of the input image. */
  typedef typename RegionType::SizeType SizeType;

  /** Pointer Type for input image. */
  typedef typename InputImageType::ConstPointer InputImagePointer;

  /** Pointer Type for the output image. */
  typedef typename OutputImageType::Pointer OutputImagePointer;
  
  /** Boundary condition type for the neighborhood iterator */
  typedef ConstantBoundaryCondition< TInputImage > ConstBoundaryConditionType;
  
  /** Neighborhood iterator type */
  typedef NeighborhoodIterator<TInputImage, ConstBoundaryConditionType> NeighborhoodIteratorType;
  
  /** Neighborhood type */
  typedef typename NeighborhoodIteratorType::NeighborhoodType NeighborhoodType;

  /** Get Skelenton by thinning image. */
  OutputImageType * GetThinning(void);

  EndPointListType & GetEndPoints(void)
  { return m_EndPoints; };

  /** ImageDimension enumeration   */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension );
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimensionCheck,
    (Concept::SameDimension<InputImageDimension, 3>));
  itkConceptMacro(SameTypeCheck,
    (Concept::SameType<InputImagePixelType, OutputImagePixelType>));
  itkConceptMacro(InputAdditiveOperatorsCheck,
    (Concept::AdditiveOperators<InputImagePixelType>));
  itkConceptMacro(InputConvertibleToIntCheck,
    (Concept::Convertible<InputImagePixelType, int>));
  itkConceptMacro(IntConvertibleToInputCheck,
    (Concept::Convertible<int, InputImagePixelType>));
  itkConceptMacro(InputIntComparableCheck,
    (Concept::Comparable<InputImagePixelType, int>));
  /** End concept checking */
#endif

protected:
  BinaryThinningImageFilter3D();
  virtual ~BinaryThinningImageFilter3D() {};

  void PrintSelf(std::ostream& os, Indent indent) const override;

  /** Compute thinning Image. */
  void GenerateData() override;

  /** Prepare data. */
  void PrepareData();

  /**  Compute thinning Image. */
  void ComputeThinImage();
  
  /**  isEulerInvariant [Lee94] */
  bool isEulerInvariant(NeighborhoodType neighbors, int *LUT);
  void fillEulerLUT(int *LUT);  
  /**  isSimplePoint [Lee94] */
  bool isSimplePoint(NeighborhoodType neighbors);
  /**  Octree_labeling [Lee94] */
  void Octree_labeling(int octant, int label, int *cube);


private:   
  BinaryThinningImageFilter3D(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  EndPointListType m_EndPoints;

}; // end of BinaryThinningImageFilter3D class

} // end namespace tube

} //end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeBinaryThinningImageFilter3D.hxx"
#endif

#endif
