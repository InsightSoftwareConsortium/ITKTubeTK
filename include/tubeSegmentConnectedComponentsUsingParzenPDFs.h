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

template< class TInputPixel, unsigned int TDimension,
  class TLabelMapPixel = unsigned char >
class SegmentConnectedComponentsUsingParzenPDFs:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef SegmentConnectedComponentsUsingParzenPDFs< TInputPixel, TDimension,
            TLabelMapPixel>                         Self;
  typedef itk::ProcessObject                        Superclass;
  typedef itk::SmartPointer< Self >                 Pointer;
  typedef itk::SmartPointer< const Self >           ConstPointer;

  typedef TInputPixel                                   InputImagePixelType;
  typedef itk::Image< InputImagePixelType, TDimension > InputImageType;
  typedef TLabelMapPixel                                LabelMapPixelType;
  typedef itk::Image< LabelMapPixelType, TDimension >   LabelMapType;
  typedef LabelMapType                                  OutputImageType;

  typedef itk::tube::PDFSegmenterParzen< InputImageType,
    LabelMapType >                                  FilterType;

  typedef typename FilterType::PDFImageType         PDFImageType;
  typedef typename FilterType::ProbabilityImageType ProbabilityImageType;
  typedef typename FilterType::VectorDoubleType     VectorDoubleType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( SegmentConnectedComponentsUsingParzenPDFs, ProcessObject );

  /** Set/Get input image */
  tubeWrapSetConstObjectMacro( Input, InputImageType, Filter );
  tubeWrapGetConstObjectMacro( Input, InputImageType, Filter );

  void SetFeatureImage( InputImageType * img );
  void AddFeatureImage( InputImageType * img );

  /** Set/Get mask image */
  tubeWrapSetObjectMacro( LabelMap, LabelMapType, Filter );
  tubeWrapGetObjectMacro( LabelMap, LabelMapType, Filter );

  tubeWrapForceSetMacro( ObjectId, LabelMapPixelType, Filter );
  tubeWrapGetMacro( ObjectId, LabelMapPixelType, Filter );
  tubeWrapAddMacro( ObjectId, LabelMapPixelType, Filter );

  tubeWrapSetMacro( VoidId, LabelMapPixelType, Filter );
  tubeWrapGetMacro( VoidId, LabelMapPixelType, Filter );

  tubeWrapSetMacro( ErodeDilateRadius, unsigned int, Filter );
  tubeWrapGetMacro( ErodeDilateRadius, unsigned int, Filter );

  tubeWrapSetMacro( DilateFirst, bool, Filter );
  tubeWrapGetMacro( DilateFirst, bool, Filter );

  tubeWrapSetMacro( HoleFillIterations, unsigned int, Filter );
  tubeWrapGetMacro( HoleFillIterations, unsigned int, Filter );

  tubeWrapSetNthMacro( ObjectPDFWeight, double, Filter );
  tubeWrapGetNthMacro( ObjectPDFWeight, double, Filter );

  tubeWrapSetMacro( ProbabilityImageSmoothingStandardDeviation, double, Filter );
  tubeWrapGetMacro( ProbabilityImageSmoothingStandardDeviation, double, Filter );

  tubeWrapSetMacro( HistogramSmoothingStandardDeviation, double, Filter );
  tubeWrapGetMacro( HistogramSmoothingStandardDeviation, double, Filter );

  tubeWrapSetMacro( ReclassifyNotObjectLabels, bool, Filter );
  tubeWrapGetMacro( ReclassifyNotObjectLabels, bool, Filter );

  tubeWrapSetMacro( ReclassifyObjectLabels, bool, Filter );
  tubeWrapGetMacro( ReclassifyObjectLabels, bool, Filter );

  tubeWrapSetMacro( ForceClassification, bool, Filter );
  tubeWrapGetMacro( ForceClassification, bool, Filter );

  tubeWrapSetMacro( BinMin, VectorDoubleType, Filter );
  tubeWrapGetMacro( BinMin, VectorDoubleType, Filter );

  tubeWrapSetMacro( BinSize, VectorDoubleType, Filter );
  tubeWrapGetMacro( BinSize, VectorDoubleType, Filter );

  tubeWrapSetNthObjectMacro( ClassPDFImage, PDFImageType, Filter );
  tubeWrapGetNthObjectMacro( ClassPDFImage, PDFImageType, Filter );

  tubeWrapCallMacro( ClassifyImages, Filter );

  tubeWrapGetMacro( NumberOfClasses, unsigned int, Filter );

  tubeWrapGetNthMacro( ClassProbabilityImage, ProbabilityImageType, Filter );

  void Update( void ) override;

  /** Get output segmentation mask */
  tubeWrapGetObjectMacro( Output, OutputImageType, Filter );

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
