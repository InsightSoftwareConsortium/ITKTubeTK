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
#ifndef tubeEnhanceTubesUsingDiscriminantAnalysis_h
#define tubeEnhanceTubesUsingDiscriminantAnalysis_h

// ITK Includes
#include "itkProcessObject.h"

// TubeTK Includes
#include "tubeWrappingMacros.h"

#include "itktubeRidgeSeedFilter.h"
#include "itktubeRidgeSeedFilterIO.h"

namespace tube
{
/** \class SegmentTubes
 *
 *  \ingroup TubeTK
 */

template< class TImage,
  class TLabelMap = itk::Image< unsigned short, TImage::ImageDimension> >
class EnhanceTubesUsingDiscriminantAnalysis:
  public itk::ProcessObject
{
public:
  /** Standard class type alias. */
  using Self = EnhanceTubesUsingDiscriminantAnalysis;
  using Superclass = itk::ProcessObject;
  using Pointer = itk::SmartPointer< Self >;
  using ConstPointer = itk::SmartPointer< const Self >;

  using ImageType = TImage;
  using PixelType = typename ImageType::PixelType;
  using IndexType = typename ImageType::IndexType;

  using LabelMapType = TLabelMap;

  using FilterType = itk::tube::RidgeSeedFilter< ImageType, LabelMapType >;
  using FilterIOType = itk::tube::RidgeSeedFilterIO< ImageType, LabelMapType >;

  using OutputImageType = typename FilterType::ProbabilityImageType;

  using ObjectIdType = typename FilterType::ObjectIdType;
  using RidgeScalesType = typename FilterType::RidgeScalesType;
  using WhitenMeansType = typename FilterType::WhitenMeansType;
  using WhitenStdDevsType = typename FilterType::WhitenStdDevsType;
  using VectorType = typename FilterType::VectorType;
  using MatrixType = typename FilterType::MatrixType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( EnhanceTubesUsingDiscriminantAnalysis, ProcessObject );

  /***/
  /***/
  /***/

  /** Set the source image. */
  tubeWrapSetConstObjectMacro( Input, ImageType, Filter );
  tubeWrapAddConstObjectMacro( Input, ImageType, Filter );
  tubeWrapGetNthConstObjectMacro( Input, ImageType, Filter );

  /** Set the labelmap image.  Is the output image after computations. */
  tubeWrapSetObjectMacro( LabelMap, LabelMapType, Filter );
  tubeWrapGetObjectMacro( LabelMap, LabelMapType, Filter );

  /***/
  /***/
  /***/

  tubeWrapSetMacro( RidgeId, ObjectIdType, Filter );
  tubeWrapGetMacro( RidgeId, ObjectIdType, Filter );

  tubeWrapSetMacro( BackgroundId, ObjectIdType, Filter );
  tubeWrapGetMacro( BackgroundId, ObjectIdType, Filter );

  tubeWrapSetMacro( UnknownId, ObjectIdType, Filter );
  tubeWrapGetMacro( UnknownId, ObjectIdType, Filter );

  tubeWrapSetMacro( IgnoreId, ObjectIdType, Filter );
  tubeWrapGetMacro( IgnoreId, ObjectIdType, Filter );

  tubeWrapSetMacro( Scales, RidgeScalesType, Filter );
  tubeWrapGetMacro( Scales, RidgeScalesType, Filter );

  tubeWrapSetMacro( InputWhitenMeans, WhitenMeansType, Filter );
  tubeWrapGetMacro( InputWhitenMeans, WhitenMeansType, Filter );

  tubeWrapSetMacro( InputWhitenStdDevs, WhitenStdDevsType, Filter );
  tubeWrapGetMacro( InputWhitenStdDevs, WhitenStdDevsType, Filter );

  tubeWrapSetMacro( OutputWhitenMeans, WhitenMeansType, Filter );
  tubeWrapGetMacro( OutputWhitenMeans, WhitenMeansType, Filter );

  tubeWrapSetMacro( OutputWhitenStdDevs, WhitenStdDevsType, Filter );
  tubeWrapGetMacro( OutputWhitenStdDevs, WhitenStdDevsType, Filter );

  tubeWrapGetMacro( NumberOfBasis, unsigned int, Filter );

  tubeWrapGetNthMacro( BasisValue, double, Filter );
  tubeWrapSetNthMacro( BasisValue, double, Filter );

  tubeWrapGetNthMacro( BasisVector, VectorType, Filter );
  tubeWrapSetNthMacro( BasisVector, VectorType, Filter );

  tubeWrapGetMacro( BasisMatrix, MatrixType, Filter );
  tubeWrapSetMacro( BasisMatrix, MatrixType, Filter );

  tubeWrapGetMacro( BasisValues, VectorType, Filter );
  tubeWrapSetMacro( BasisValues, VectorType, Filter );

  tubeWrapGetNthObjectMacro( ClassProbabilityImage, OutputImageType, Filter );
  tubeWrapGetNthObjectMacro( ClassLikelihoodRatioImage, OutputImageType, Filter );

  /** Load parameters of tube extraction from a file */
  void LoadParameterFile( const std::string & filename )
  { FilterIOType teReader( m_Filter );
    teReader.Read( filename.c_str() ); };

  void SaveParameterFile( const std::string & filename )
  { FilterIOType teWriter( m_Filter );
    teWriter.Write( filename.c_str() ); };

  /***/
  /***/
  /***/

  tubeWrapSetMacro( SeedTolerance, double, Filter );
  tubeWrapGetMacro( SeedTolerance, double, Filter );

  tubeWrapSetMacro( Skeletonize, bool, Filter );
  tubeWrapGetMacro( Skeletonize, bool, Filter );

  tubeWrapSetMacro( UseIntensityOnly, bool, Filter );
  tubeWrapGetMacro( UseIntensityOnly, bool, Filter );

  tubeWrapSetMacro( UseFeatureMath, bool, Filter );
  tubeWrapGetMacro( UseFeatureMath, bool, Filter );

  tubeWrapSetMacro( TrainClassifier, bool, Filter );
  tubeWrapGetMacro( TrainClassifier, bool, Filter );

  tubeWrapUpdateMacro( Filter );

  tubeWrapCallMacro( ClassifyImages, Filter );

  tubeWrapGetObjectMacro( Output, LabelMapType, Filter);
  tubeWrapGetObjectMacro( OutputSeedScales, OutputImageType, Filter);

protected:
  EnhanceTubesUsingDiscriminantAnalysis( void );
  ~EnhanceTubesUsingDiscriminantAnalysis() {};
  void PrintSelf( std::ostream & os, itk::Indent indent ) const override;

private:
  EnhanceTubesUsingDiscriminantAnalysis( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) 
    override {};

  typename FilterType::Pointer           m_Filter;

};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeEnhanceTubesUsingDiscriminantAnalysis.hxx"
#endif

#endif // End !defined( __EnhanceTubesUsingDiscriminantAnalysis_h )
