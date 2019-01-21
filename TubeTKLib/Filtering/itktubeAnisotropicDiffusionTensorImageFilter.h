/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

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

#ifndef __itktubeAnisotropicDiffusionTensorImageFilter_h
#define __itktubeAnisotropicDiffusionTensorImageFilter_h

#include "itktubeAnisotropicDiffusionTensorFunction.h"

#include <itkFiniteDifferenceImageFilter.h>

namespace itk
{

namespace tube
{

/** \class AnisotropicDiffusionTensorImageFilter
 * \brief This is a superclass for filters that iteratively enhance edges in
 *        an image by solving a non-linear diffusion equation.
 *
 * \warning Does not handle image directions.  Re-orient images to axial
 * ( direction cosines = identity matrix ) before using this function.
 *
 * \sa AnisotropicEdgeEnhancementDiffusionImageFilter
 * \sa AnisotropicCoherenceEnhancingDiffusionImageFilter
 * \sa AnisotropicHybridDiffusionImageFilter
 *
 * \ingroup FiniteDifferenceFunctions
 * \ingroup Functions
 */
template< class TInputImage, class TOutputImage >
class AnisotropicDiffusionTensorImageFilter
  : public FiniteDifferenceImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs */
  typedef AnisotropicDiffusionTensorImageFilter                   Self;
  typedef FiniteDifferenceImageFilter<TInputImage, TOutputImage>  Superclass;
  typedef SmartPointer< Self >                                    Pointer;
  typedef SmartPointer< const Self >                              ConstPointer;

  // itkNewMacro( Self );  // Not included since pure virtual

  /** Run-time type information ( and related methods ) */
  itkTypeMacro( AnisotropicDiffusionTensorImageFilter,
    FiniteDifferenceImageFiler );

  /** Convenient typedefs */
  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename Superclass::PixelType       PixelType;

  /** Dimensionality of input and output data is assumed to be the same.
   * It is inherited from the superclass. */
  itkStaticConstMacro( ImageDimension, unsigned int,
    Superclass::ImageDimension );

  /** Type of associated function, with associated typedefs */
  typedef AnisotropicDiffusionTensorFunction< InputImageType >
      FiniteDifferenceFunctionType;
  typedef typename FiniteDifferenceFunctionType::DiffusionTensorType
      TensorPixelType;
  typedef typename FiniteDifferenceFunctionType::DiffusionTensorImageType
      DiffusionTensorImageType;

  // Define the type for storing the eigenvalues
  typedef FixedArray< double, ImageDimension >   EigenValueArrayType;

  // Declare the types of the output images
  typedef Image< EigenValueArrayType, ImageDimension >
      EigenAnalysisOutputImageType;

  /** The value type of a time step.  Inherited from the superclass. */
  typedef typename Superclass::TimeStepType TimeStepType;

  /** The container type for the update buffer. */
  typedef OutputImageType UpdateBufferType;

  /** Define diffusion image neighborhood type */
  typedef typename
    FiniteDifferenceFunctionType::DiffusionTensorNeighborhoodType
      DiffusionTensorNeighborhoodType;

  /** Set/Get Macro for diffusion tensor image filter parameters */
  itkSetMacro( TimeStep, double );
  itkGetMacro( TimeStep, double );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( OutputTimesDoubleCheck,
    ( Concept::MultiplyOperator< PixelType, double > ) );
  itkConceptMacro( OutputAdditiveOperatorsCheck,
    ( Concept::AdditiveOperators< PixelType > ) );
  itkConceptMacro( InputConvertibleToOutputCheck,
    ( Concept::Convertible< typename TInputImage::PixelType, PixelType > ) );
  /** End concept checking */
#endif

protected:
  AnisotropicDiffusionTensorImageFilter( void );
 ~AnisotropicDiffusionTensorImageFilter( void ) {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

  /* overloaded GenerateData method */
  virtual void GenerateData( void );

  /** A simple method to copy the data from the input to the output. ( Supports
   * "read-only" image adaptors in the case where the input image type converts
   * to a different output image type. )  */
  virtual void CopyInputToOutput( void );

  /** This method applies changes from the m_UpdateBuffer to the output using
   * the ThreadedApplyUpdate() method and a multithreading mechanism.  "dt" is
   * the time step to use for the update of each pixel. */
  virtual void ApplyUpdate( const TimeStepType& dt );

  /** Method to allow subclasses to get direct access to the update
   * buffer */
  virtual UpdateBufferType* GetUpdateBuffer( void )
    { return m_UpdateBuffer; }

  /** This method populates an update buffer with changes for each pixel in the
   * output using the ThreadedCalculateChange() method and a multithreading
   * mechanism. Returns value is a time step to be used for the update. */
  virtual TimeStepType CalculateChange( void );

  /** This method allocates storage in m_UpdateBuffer.  It is called from
   * Superclass::GenerateData(). */
  virtual void AllocateUpdateBuffer( void );

  /** This method allocates storage for the diffusion tensor image */
  void AllocateDiffusionTensorImage( void );

  /** Update diffusion tensor image */
  void virtual UpdateDiffusionTensorImage( void ) = 0;

  /** The type of region used for multithreading */
  typedef typename UpdateBufferType::RegionType ThreadRegionType;

  /** The type of region used for multithreading */
  typedef typename DiffusionTensorImageType::RegionType
      ThreadDiffusionTensorImageRegionType;

  typedef typename DiffusionTensorImageType::Pointer
      DiffusionTensorImagePointerType;

  /**  Does the actual work of updating the output from the UpdateContainer
   *   over an output region supplied by the multithreading mechanism.
   *  \sa ApplyUpdate
   *  \sa ApplyUpdateThreaderCallback */
  virtual void ThreadedApplyUpdate( TimeStepType dt,
    const ThreadRegionType &regionToProcess,
    const ThreadDiffusionTensorImageRegionType &diffusionRegionToProcess,
    ThreadIdType threadId );

  /** Does the actual work of calculating change over a region supplied by
   * the multithreading mechanism.
   * \sa CalculateChange
   * \sa CalculateChangeThreaderCallback */
  virtual TimeStepType ThreadedCalculateChange(
    const ThreadRegionType &regionToProcess,
    const ThreadDiffusionTensorImageRegionType &diffusionRegionToProcess,
    ThreadIdType threadId );

  /** Prepare for the iteration process. */
  virtual void InitializeIteration( void );

  DiffusionTensorImagePointerType GetDiffusionTensorImage( void );

private:
  //purposely not implemented
  AnisotropicDiffusionTensorImageFilter( const Self& );
  void operator=( const Self& ); //purposely not implemented

  /** Structure for passing information into static callback methods.  Used in
   * the subclasses' threading mechanisms. */
  struct DenseFDThreadStruct
    {
    AnisotropicDiffusionTensorImageFilter *Filter;
    TimeStepType TimeStep;
    std::vector< TimeStepType > TimeStepList;
    std::vector< bool > ValidTimeStepList;

    }; // End struct DenseFDThreadStruct

  /** This callback method uses ImageSource::SplitRequestedRegion to acquire an
   * output region that it passes to ThreadedApplyUpdate for processing. */
  static ITK_THREAD_RETURN_TYPE ApplyUpdateThreaderCallback( void *arg );

  /** This callback method uses SplitUpdateContainer to acquire a region
   * which it then passes to ThreadedCalculateChange for processing. */
  static ITK_THREAD_RETURN_TYPE CalculateChangeThreaderCallback( void *arg );

  typename DiffusionTensorImageType::Pointer            m_DiffusionTensorImage;

  /** The buffer that holds the updates for an iteration of the algorithm. */
  typename UpdateBufferType::Pointer                    m_UpdateBuffer;

  TimeStepType                                          m_TimeStep;

}; // End class AnisotropicDiffusionTensorImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeAnisotropicDiffusionTensorImageFilter.hxx"
#endif

#endif // End !defined( __itktubeAnisotropicDiffusionTensorImageFilter_h )
