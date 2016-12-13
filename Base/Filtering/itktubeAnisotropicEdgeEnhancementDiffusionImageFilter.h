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

#ifndef __itktubeAnisotropicEdgeEnhancementDiffusionImageFilter_h
#define __itktubeAnisotropicEdgeEnhancementDiffusionImageFilter_h

#include "itktubeAnisotropicDiffusionTensorFunction.h"
#include "itktubeAnisotropicDiffusionTensorImageFilter.h"
#include "itktubeStructureTensorRecursiveGaussianImageFilter.h"
#include "itktubeSymmetricEigenVectorAnalysisImageFilter.h"

#include <itkDiffusionTensor3D.h>
#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkMultiThreader.h>
#include <itkSymmetricEigenAnalysisImageFilter.h>

namespace itk
{

namespace tube
{
/** \class AnisotropicEdgeEnhancementDiffusionImageFilter
 *  This class is an implementation of Edge-enhancing diffusion:
 *   Mendrik et al., Noise reduction in computed tomography scans using 3-D
 *   anisotropic hybrid diffusion with continuous switch. IEEE Transactions on
 *   Medical Imaging 28( 10 ), pp. 1585-1594, 2009.
 *
 * \warning Does not handle image directions.  Re-orient images to axial
 * ( direction cosines = identity matrix ) before using this function.
 *
 * \sa AnisotropicDiffusionTensorImageFilter
 * \sa AnisotropicCoherenceEnhancingDiffusionImageFilter
 * \sa AnisotropicHybridDiffusionImageFilter
 *
 * \ingroup FiniteDifferenceFunctions
 * \ingroup Functions
 */


template< class TInputImage, class TOutputImage >
class AnisotropicEdgeEnhancementDiffusionImageFilter
  : public AnisotropicDiffusionTensorImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs */
  typedef AnisotropicEdgeEnhancementDiffusionImageFilter Self;

  typedef AnisotropicDiffusionTensorImageFilter<TInputImage, TOutputImage>
                                                           Superclass;

  typedef SmartPointer< Self >                             Pointer;
  typedef SmartPointer< const Self >                       ConstPointer;


  /** Method for creation through the object factory */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ) */
  itkTypeMacro( AnisotropicEdgeEnhancementDiffusionImageFilter,
                ImageToImageFilter );

  /** Convenient typedefs */
  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename Superclass::PixelType       PixelType;

  typedef typename Superclass::DiffusionTensorImageType
                                                DiffusionTensorImageType;

  // Structure tensor type
  typedef StructureTensorRecursiveGaussianImageFilter < InputImageType >
                                                StructureTensorFilterType;

  /** Dimensionality of input and output data is assumed to be the same.
   * It is inherited from the superclass. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       Superclass::ImageDimension );

  typedef itk::Matrix<double, ImageDimension, ImageDimension> MatrixType;

  // Define image of matrix pixel type
  typedef itk::Image< MatrixType, ImageDimension>  OutputMatrixImageType;

  // Define the symmetric tensor pixel type
  typedef itk::SymmetricSecondRankTensor< double, ImageDimension>
                                                         TensorPixelType;
  typedef itk::Image< TensorPixelType, ImageDimension>
                                                         TensorImageType;

   // Define the type for storing the eigenvalue
  typedef itk::FixedArray< double, ImageDimension >      EigenValueArrayType;

  // Declare the types of the output images
  typedef itk::Image< EigenValueArrayType, ImageDimension >
    EigenAnalysisOutputImageType;

  /** The container type for the update buffer. */
  typedef OutputImageType UpdateBufferType;

  /** Define diffusion image nbd type */
  typedef typename Superclass::DiffusionTensorNeighborhoodType
    DiffusionTensorNeighborhoodType;

  /** Set the contrast parameter */
  itkSetMacro( ContrastParameterLambdaE, double );

  /** Set threshold parameter C */
  itkSetMacro( ThresholdParameterC, double );

  /** Set the sigma and outer sigma value for structure tensor computation */
  itkSetMacro( Sigma, double );
  itkSetMacro( SigmaOuter, double );

  /**Get methods */
  itkGetMacro( ContrastParameterLambdaE, double );
  itkGetMacro( ThresholdParameterC, double );
  itkGetMacro( Sigma, double );
  itkGetMacro( SigmaOuter, double );


protected:
  AnisotropicEdgeEnhancementDiffusionImageFilter( void );
 ~AnisotropicEdgeEnhancementDiffusionImageFilter( void ) {}
  void PrintSelf( std::ostream& os, Indent indent ) const;

  /** Update diffusion tensor image */
  void virtual UpdateDiffusionTensorImage( void );

private:
  //purposely not implemented
  AnisotropicEdgeEnhancementDiffusionImageFilter( const Self& );
  void operator=( const Self& ); //purposely not implemented

  double    m_ContrastParameterLambdaE;
  double    m_ThresholdParameterC;
  double    m_Sigma;
  double    m_SigmaOuter;

}; // End class AnisotropicEdgeEnhancementDiffusionImageFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeAnisotropicEdgeEnhancementDiffusionImageFilter.hxx"
#endif

#endif
// End !defined( __itktubeAnisotropicEdgeEnhancementDiffusionImageFilter_h )
