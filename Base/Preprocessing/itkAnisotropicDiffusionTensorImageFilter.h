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
#ifndef __itkAnisotropicDiffusionTensorImageFilter_h
#define __itkAnisotropicDiffusionTensorImageFilter_h

#include "itkFiniteDifferenceImageFilter.h"
#include "itkAnisotropicDiffusionTensorFunction.h"
#include "itkMultiThreader.h"
#include "itkDiffusionTensor3D.h"
#include "itkSymmetricEigenAnalysisImageFilter.h"
#include "itkSymmetricEigenVectorAnalysisImageFilter.h"

namespace itk {
/** \class AnisotropicDiffusionTensorImageFilter
 * \brief This is a superclass for filters that iteratively enhance edge in 
 *        an image by solving non-linear diffusion equation.
 *
 * 
 * \sa AnisotropicEdgeEnhancementDiffusionImageFilter
 * \sa AnisotropicCoherenceEnhancingDiffusionImageFilter
 *
 * \ingroup FiniteDifferenceFunctions
 * \ingroup Functions
 */


template <class TInputImage, class TOutputImage>
class ITK_EXPORT AnisotropicDiffusionTensorImageFilter  
  : public FiniteDifferenceImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs */
  typedef AnisotropicDiffusionTensorImageFilter Self;

  typedef FiniteDifferenceImageFilter<TInputImage, TOutputImage> 
                                                           Superclass;

  typedef SmartPointer<Self>                               Pointer;
  typedef SmartPointer<const Self>                         ConstPointer;
 

  /** Run-time type information (and related methods) */
  itkTypeMacro(AnisotropicDiffusionTensorImageFilter,
                                                ImageToImageFilter );
  
  /** Convenient typedefs */
  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename Superclass::PixelType       PixelType;

  typedef itk::Image< DiffusionTensor3D< double > , 3 > 
                                                DiffusionTensorImageType;

  /** Dimensionality of input and output data is assumed to be the same.
   * It is inherited from the superclass. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  typedef AnisotropicDiffusionTensorFunction<InputImageType>  
                                                  FiniteDifferenceFunctionType;
  
  typedef itk::Image< double, 3 >               VesselnessOutputImageType;

  typedef itk::Matrix<double, ImageDimension, ImageDimension> MatrixType;

  // Define image of matrix pixel type 
  typedef itk::Image< MatrixType, ImageDimension>  OutputMatrixImageType;

  // Define the symmetric tensor pixel type
  typedef itk::SymmetricSecondRankTensor< double, ImageDimension> 
                                                         TensorPixelType;
  typedef itk::Image< TensorPixelType, ImageDimension>  
                                                         TensorImageType;

   // Define the type for storing the eigen-value
  typedef itk::FixedArray< double, ImageDimension >      EigenValueArrayType;
  
  // Declare the types of the output images
  typedef itk::Image< EigenValueArrayType, ImageDimension >  
                                                  EigenAnalysisOutputImageType;
  
  /** The value type of a time step.  Inherited from the superclass. */
  typedef typename Superclass::TimeStepType TimeStepType;

  /** The container type for the update buffer. */
  typedef OutputImageType UpdateBufferType;

  /** Define diffusion image nbd type */
  typedef typename FiniteDifferenceFunctionType::DiffusionTensorNeighborhoodType
                                               DiffusionTensorNeighborhoodType;


  /** Set/Get Macro for VED parameters */
  itkSetMacro( TimeStep, double ); 

  itkGetMacro( TimeStep, double ); 

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(OutputTimesDoubleCheck,
    (Concept::MultiplyOperator<PixelType, double>));
  itkConceptMacro(OutputAdditiveOperatorsCheck,
    (Concept::AdditiveOperators<PixelType>));
  itkConceptMacro(InputConvertibleToOutputCheck,
    (Concept::Convertible<typename TInputImage::PixelType, PixelType>));
  /** End concept checking */
#endif

protected:
  AnisotropicDiffusionTensorImageFilter();
 ~AnisotropicDiffusionTensorImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  /* overloaded GenerateData method */
  virtual void GenerateData(); 

  /** A simple method to copy the data from the input to the output. ( Supports
   * "read-only" image adaptors in the case where the input image type converts
   * to a different output image type. )  */
  virtual void CopyInputToOutput();

  /** This method applies changes from the m_UpdateBuffer to the output using
   * the ThreadedApplyUpdate() method and a multithreading mechanism.  "dt" is
   * the time step to use for the update of each pixel. */
  virtual void ApplyUpdate(TimeStepType dt);

  /** Method to allow subclasses to get direct access to the update
   * buffer */
  virtual UpdateBufferType* GetUpdateBuffer()
    { return m_UpdateBuffer; }

  /** This method populates an update buffer with changes for each pixel in the
   * output using the ThreadedCalculateChange() method and a multithreading
   * mechanism. Returns value is a time step to be used for the update. */
  virtual TimeStepType CalculateChange();

  /** This method allocates storage in m_UpdateBuffer.  It is called from
   * Superclass::GenerateData(). */
  virtual void AllocateUpdateBuffer();

  /** This method allocates storage for the diffusion tensor image */
  void AllocateDiffusionTensorImage();
 
  /** Update diffusion tensor image */
  void virtual UpdateDiffusionTensorImage() = 0;
 
  /** The type of region used for multithreading */
  typedef typename UpdateBufferType::RegionType ThreadRegionType;

  /** The type of region used for multithreading */
  typedef typename DiffusionTensorImageType::RegionType 
                                        ThreadDiffusionTensorImageRegionType;

  typedef typename DiffusionTensorImageType::Pointer DiffusionTensorImagePointerType;

  /**  Does the actual work of updating the output from the UpdateContainer 
   *   over an output region supplied by the multithreading mechanism.
   *  \sa ApplyUpdate
   *  \sa ApplyUpdateThreaderCallback */ 
  virtual
  void ThreadedApplyUpdate(
                TimeStepType dt,
                const ThreadRegionType &regionToProcess,
                const ThreadDiffusionTensorImageRegionType &diffusionRegionToProcess,
                int threadId);

  /** Does the actual work of calculating change over a region supplied by
   * the multithreading mechanism.
   * \sa CalculateChange
   * \sa CalculateChangeThreaderCallback */
  virtual
  TimeStepType ThreadedCalculateChange(
               const ThreadRegionType &regionToProcess,
               const ThreadDiffusionTensorImageRegionType &diffusionRegionToProcess,
               int threadId);

  /** Prepare for the iteration process. */
  virtual void InitializeIteration();

  DiffusionTensorImagePointerType GetDiffusionTensorImage();

private:
  //purposely not implemented
  AnisotropicDiffusionTensorImageFilter(const Self&); 
  void operator=(const Self&); //purposely not implemented

  /** Structure for passing information into static callback methods.  Used in
   * the subclasses' threading mechanisms. */
  struct DenseFDThreadStruct
    {
    AnisotropicDiffusionTensorImageFilter *Filter;
    TimeStepType TimeStep;
    TimeStepType *TimeStepList;
    bool *ValidTimeStepList;
    };
    
  /** This callback method uses ImageSource::SplitRequestedRegion to acquire an
   * output region that it passes to ThreadedApplyUpdate for processing. */
  static ITK_THREAD_RETURN_TYPE ApplyUpdateThreaderCallback( void *arg );
  
  /** This callback method uses SplitUpdateContainer to acquire a region
   * which it then passes to ThreadedCalculateChange for processing. */
  static ITK_THREAD_RETURN_TYPE CalculateChangeThreaderCallback( void *arg );
 
  typename DiffusionTensorImageType::Pointer            m_DiffusionTensorImage;

  /** The buffer that holds the updates for an iteration of the algorithm. */
  typename UpdateBufferType::Pointer m_UpdateBuffer;

  TimeStepType                                          m_TimeStep;

};
  

}// end namespace itk

#if ITK_TEMPLATE_TXX
# include "itkAnisotropicDiffusionTensorImageFilter.txx"
#endif

#endif
