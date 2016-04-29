/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/
#ifndef __tubeConvertSpatialGraphToImage_h
#define __tubeConvertSpatialGraphToImage_h

#include "itktubeConvertSpatialGraphToImageFilter.h"
#include "itkObject.h"


namespace tube
{
/** \class ConvertSpatialGraphToImage
 *
 *  \ingroup TubeTKITK
 */

template< typename TInputImage, typename TOutputImage >
class ConvertSpatialGraphToImage:
  public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef ConvertSpatialGraphToImage      Self;
  typedef itk::SmartPointer< Self >       Pointer;
  typedef itk::SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ConvertSpatialGraphToImage, Object);

  typename TOutputImage::Pointer GetAdjacencyMatrixImage();

  typename TOutputImage::Pointer GetBranchnessImage();

  typename TOutputImage::Pointer GetRadiusImage();

  typename TOutputImage::Pointer GetCentralityImage();

  void SetAdjacencyMatrix( vnl_matrix< double > );

  void SetBranchnessVector( vnl_vector< double > );

  void SetRadiusVector( vnl_vector< double > );

  void SetCentralityVector( vnl_vector< double > );

  void SetInput( const TInputImage *inputImage );

  void Update();

  typename TOutputImage::Pointer GetOutput();

protected:
  ConvertSpatialGraphToImage( void );
  ~ConvertSpatialGraphToImage() {}
  void PrintSelf(std::ostream & os, itk::Indent indent) const;

private:
  /** itkConvertSpatialGraphToImageFilter parameters **/
  ConvertSpatialGraphToImage(const Self &);
  void operator=(const Self &);

  typedef itk::tube::ConvertSpatialGraphToImageFilter
    < TInputImage, TOutputImage > ConvertSpatialGraphToImageFilterType;
  typename ConvertSpatialGraphToImageFilterType::Pointer
    m_ConvertSpatialGraphToImageFilter;

};
} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeConvertSpatialGraphToImage.hxx"
#endif

#endif // End !defined( __tubeConvertSpatialGraphToImage_h )
