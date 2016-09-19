/*=========================================================================

Library:   TubeTK/VTree3D

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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

#ifndef __itktubeComputeSegmentTubesParameters_h
#define __itktubeComputeSegmentTubesParameters_h

#include "itkImage.h"
#include "itkObject.h"
#include <algorithm>
#include <vector>
#include <itktubeRidgeExtractor.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkContinuousIndex.h>

namespace itk
{

namespace tube
{

/**
 * This class computes the parameters for the segmentation
 *
 * \sa ComputeSegmentTubesParameters
 */

template< class TPixel, unsigned int VDimension >
class ComputeSegmentTubesParameters : public Object
{
public:

  /**
   * Standard self typedef */
  typedef ComputeSegmentTubesParameters Self;
  typedef Object                        Superclass;
  typedef SmartPointer< Self >          Pointer;
  typedef SmartPointer< const Self >    ConstPointer;

  typedef Image< TPixel, VDimension >         InputImageType;
  typedef typename InputImageType::Pointer    InputImagePointer;
  typedef int                                 MaskPixelType;
  typedef Image< MaskPixelType, VDimension >  MaskImageType;
  typedef typename MaskImageType::Pointer     MaskImagePointer;
  typedef float                               ScalePixelType;
  typedef Image< ScalePixelType, VDimension > ScaleImageType;
  typedef typename ScaleImageType::Pointer    ScaleImagePointer;

  typedef vnl_vector< double >                MetricVectorType;
  typedef std::vector< MetricVectorType >     SampleListType;

  typedef itk::ContinuousIndex< double, VDimension > IndexType;
  typedef std::vector< IndexType >                   IndexListType;

  typedef itk::tube::RidgeExtractor< InputImageType >
  RidgeExtractorType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ComputeSegmentTubesParameters, Object );

  /** Set/Get the input image */
  itkSetObjectMacro( InputImage, InputImageType );
  itkGetObjectMacro( InputImage, InputImageType );

  /** Set/Get the mask input image */
  itkSetObjectMacro( MaskInputImage, MaskImageType );
  itkGetObjectMacro( MaskInputImage, MaskImageType );

  /** Set/Get the scale input image */
  itkSetObjectMacro( ScaleInputImage, ScaleImageType );
  itkGetObjectMacro( ScaleInputImage, ScaleImageType );

  /** Set/Get Mask BackGround Id */
  itkSetMacro( MaskBackGroundId, int );
  itkGetMacro( MaskBackGroundId, int );

  /** Set/Get Mask BackGround Id */
  itkSetMacro( MaskTubeId, int );
  itkGetMacro( MaskTubeId, int );

  /** Set Parameter File */
  itkSetMacro( ParameterFile, std::string );

  /** Get Seed Data List */
  std::vector< vnl_vector< double > >
  GetSeedData( void );

  /** Get Tube Data List */
  std::vector< vnl_vector< double > >
  GetTubeData( void );

  /** Get Bkg Data List */
  std::vector< vnl_vector< double > >
  GetBkgData( void );

  /** Get Seed Index List */
  std::vector< itk::ContinuousIndex< double, VDimension > >
  GetSeedDataIndexList( void );

  /** Get Tube Index List */
  std::vector< itk::ContinuousIndex< double, VDimension > >
  GetTubeDataIndexList( void );

  /** Get Bkg Index List */
  std::vector< itk::ContinuousIndex< double, VDimension > >
  GetBkgDataIndexList( void );

  void Update();
protected:

  ComputeSegmentTubesParameters( void );
  virtual ~ComputeSegmentTubesParameters( void );
  void PrintSelf( std::ostream & os, Indent indent ) const;

private:

  ComputeSegmentTubesParameters( const Self& );
  void operator=( const Self& );

  InputImagePointer    m_InputImage;
  MaskImagePointer     m_MaskInputImage;
  ScaleImagePointer    m_ScaleInputImage;
  SampleListType       m_SeedData;
  SampleListType       m_TubeData;
  SampleListType       m_BkgData;
  IndexListType        m_SeedDataIndexList;
  IndexListType        m_TubeDataIndexList;
  IndexListType        m_BkgDataIndexList;
  int                  m_MaskBackGroundId;
  int                  m_MaskTubeId;
  std::string          m_ParameterFile;

}; // End class SegmentTubes

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeComputeSegmentTubesParameters.hxx"
#endif

#endif // End !defined(__itktubeComputeSegmentTubesParameters_h)
