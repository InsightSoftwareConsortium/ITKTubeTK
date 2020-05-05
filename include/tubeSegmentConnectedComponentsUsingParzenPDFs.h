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
#ifndef __tubeSegmentConnectedComponentsUsingParzenPDFs_h
#define __tubeSegmentConnectedComponentsUsingParzenPDFs_h

// ITK includes
#include "itkProcessObject.h"

// TubeTK includes
#include "tubeWrappingMacros.h"

#include "itktubePDFSegmenterParzen.h"


namespace tube
{
/** \class SegmentConnectedComponentsUsingParzenPDFs
 *
 *  \ingroup TubeTK
 */

template< class TImage, class TLabelMap=itk::Image<
  typename TImage::PixelType, TImage::ImageDimension> >
class SegmentConnectedComponentsUsingParzenPDFs:
  public itk::ProcessObject
{
public:
  typedef TImage                                InputImageType;
  typedef typename TImage::PixelType            InputImagePixelType;

  typedef TLabelMap                             LabelMapType;
  typedef typename TLabelMap::PixelType         LabelMapPixelType;

  /** Standard class typedefs. */
  typedef SegmentConnectedComponentsUsingParzenPDFs  Self;
  typedef itk::ProcessObject                         Superclass;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;


  typedef itk::tube::PDFSegmenterParzen< InputImageType,
    LabelMapType >                                  FilterType;

  typedef typename FilterType::PDFImageType            PDFImageType;
  //typedef typename FilterType::LabeledFeatureSpaceType LabeledFeatureSpaceType;

  typedef typename FilterType::ProbabilityImageType  ProbabilityImageType;
  typedef typename FilterType::ProbabilityVectorType ProbabilityVectorType;
  typedef typename FilterType::FeatureVectorType     FeatureVectorType;
  typedef typename FilterType::VectorDoubleType      VectorDoubleType;
  typedef typename FilterType::VectorUIntType        VectorUIntType;
  typedef typename FilterType::ObjectIdListType      ObjectIdListType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( SegmentConnectedComponentsUsingParzenPDFs, ProcessObject );

  /** Set/Get input image */
  void SetFeatureImage( InputImageType * img );
  void AddFeatureImage( InputImageType * img );

  /** Set/Get mask image */
  tubeWrapSetObjectMacro( InputLabelMap, LabelMapType, Filter );
  tubeWrapGetObjectMacro( InputLabelMap, LabelMapType, Filter );

  tubeWrapForceSetMacro( ObjectId, LabelMapPixelType, Filter );
  tubeWrapAddMacro( ObjectId, LabelMapPixelType, Filter );
  tubeWrapGetMacro( ObjectId, ObjectIdListType, Filter );
  tubeWrapGetNthObjectMacro( ObjectId, LabelMapPixelType, Filter );

  tubeWrapSetMacro( VoidId, LabelMapPixelType, Filter );
  tubeWrapGetMacro( VoidId, LabelMapPixelType, Filter );

  tubeWrapSetMacro( ErodeDilateRadius, unsigned int, Filter );
  tubeWrapGetMacro( ErodeDilateRadius, unsigned int, Filter );

  tubeWrapSetMacro( DilateFirst, bool, Filter );
  tubeWrapGetMacro( DilateFirst, bool, Filter );

  tubeWrapSetMacro( HoleFillIterations, unsigned int, Filter );
  tubeWrapGetMacro( HoleFillIterations, unsigned int, Filter );

  tubeWrapSetNthMacro( ObjectPDFWeight, double, Filter );
  tubeWrapGetNthObjectMacro( ObjectPDFWeight, double, Filter );
  tubeWrapSetMacro( ObjectPDFWeight, VectorDoubleType, Filter );
  tubeWrapGetMacro( ObjectPDFWeight, VectorDoubleType, Filter );

  tubeWrapSetMacro( ProbabilityImageSmoothingStandardDeviation, double, Filter );
  tubeWrapGetMacro( ProbabilityImageSmoothingStandardDeviation, double, Filter );

  tubeWrapSetMacro( HistogramSmoothingStandardDeviation, double, Filter );
  tubeWrapGetMacro( HistogramSmoothingStandardDeviation, double, Filter );

  tubeWrapSetMacro( OutlierRejectPortion, double, Filter );
  tubeWrapGetMacro( OutlierRejectPortion, double, Filter );

  tubeWrapSetMacro( ReclassifyNotObjectLabels, bool, Filter );
  tubeWrapGetMacro( ReclassifyNotObjectLabels, bool, Filter );

  tubeWrapSetMacro( ReclassifyObjectLabels, bool, Filter );
  tubeWrapGetMacro( ReclassifyObjectLabels, bool, Filter );

  tubeWrapSetMacro( ForceClassification, bool, Filter );
  tubeWrapGetMacro( ForceClassification, bool, Filter );

  tubeWrapSetMacro( NumberOfBinsPerFeature, VectorUIntType, Filter );
  tubeWrapGetMacro( NumberOfBinsPerFeature, VectorUIntType, Filter );

  tubeWrapSetMacro( BinMin, VectorDoubleType, Filter );
  tubeWrapGetMacro( BinMin, VectorDoubleType, Filter );

  tubeWrapSetMacro( BinSize, VectorDoubleType, Filter );
  tubeWrapGetMacro( BinSize, VectorDoubleType, Filter );

  tubeWrapSetNthObjectMacro( ClassPDFImage, PDFImageType, Filter );
  tubeWrapGetNthObjectMacro( ClassPDFImage, PDFImageType, Filter );
  
  tubeWrapCallMacro( GenerateLabeledFeatureSpace, Filter );
  
  //tubeWrapSetObjectMacro( LabeledFeatureSpace, LabeledFeatureSpaceType, Filter );
  //tubeWrapGetObjectMacro( LabeledFeatureSpace, LabeledFeatureSpaceType, Filter );
  
  ProbabilityVectorType GetProbabilityVector( const FeatureVectorType & fv )
    const 
  { return m_Filter->GetProbabilityVector( fv ); }; 

  tubeWrapCallMacro( ClassifyImages, Filter );

  tubeWrapGetMacro( NumberOfClasses, unsigned int, Filter );

  tubeWrapGetNthObjectMacro( ClassProbabilityImage, ProbabilityImageType,
    Filter );

  void Update( void ) override;

  /** Get output segmentation mask */
  tubeWrapGetObjectMacro( OutputLabelMap, LabelMapType, Filter );

protected:
  SegmentConnectedComponentsUsingParzenPDFs( void );
  ~SegmentConnectedComponentsUsingParzenPDFs() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const override;

private:
  /** itkSegmentConnectedComponentsUsingParzenPDFs parameters **/
  SegmentConnectedComponentsUsingParzenPDFs( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) override
    {};

  typedef itk::tube::FeatureVectorGenerator< InputImageType >
                                                  FeatureVectorGeneratorType;

  typename FilterType::Pointer                    m_Filter;
  typename FeatureVectorGeneratorType::Pointer    m_FVGenerator;

};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeSegmentConnectedComponentsUsingParzenPDFs.hxx"
#endif

#endif // End !defined( __tubeSegmentConnectedComponentsUsingParzenPDFs_h )
