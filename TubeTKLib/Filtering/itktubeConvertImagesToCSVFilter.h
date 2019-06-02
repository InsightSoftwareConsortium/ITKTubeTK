/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
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
#ifndef __itktubeConvertImagesToCSVFilter_h
#define __itktubeConvertImagesToCSVFilter_h

#include <itkProcessObject.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageRegionIterator.h>
#include <itkSimpleDataObjectDecorator.h>

#include "tubeMessage.h"

namespace itk
{
namespace tube
{
/** \class ConvertImagesToCSV
 */

 template< class TInputImage, class TInputMask >
 class ConvertImagesToCSVFilter : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ConvertImagesToCSVFilter              Self;
  typedef ProcessObject                         Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;

  typedef TInputImage                                InputImageType;
  typedef typename InputImageType::Pointer           ImagePointer;
  typedef typename InputImageType::PixelType         InputPixelType;
  typedef typename InputImageType::IndexType         IndexType;
  typedef ImageFileReader< InputImageType >          ReaderType;
  typedef ImageRegionConstIterator< InputImageType > InputImageIteratorType;

  typedef TInputMask                                 InputMaskType;
  typedef typename InputMaskType::Pointer            MaskPointer;
  typedef typename InputMaskType::PixelType          MaskPixelType;
  typedef typename InputMaskType::IndexType          MaskIndexType;
  typedef ImageFileReader< InputMaskType >           MaskReaderType;
  typedef ImageRegionConstIterator< InputMaskType >  MaskIteratorType;

  typedef vnl_matrix< InputPixelType >                VnlMatrixType;
  typedef SimpleDataObjectDecorator< VnlMatrixType >  OutputType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( ConvertImagesToCSVFilter, ProcessObject );

  /** ImageDimension constants */
  itkStaticConstMacro( ImageDimension, unsigned int,
                      TInputMask::ImageDimension );

  itkSetObjectMacro( InputMask, InputMaskType );
  itkGetObjectMacro( InputMask, InputMaskType );

  OutputType* GetOutput();

  itkGetMacro( Stride, unsigned int );

  itkSetClampMacro( Stride, unsigned int, 1,
    std::numeric_limits<unsigned int>::max() );

  itkGetMacro( NumImages, unsigned int );
  itkSetMacro( NumImages, unsigned int );

  itkSetMacro( NumberRows, unsigned int );
  itkGetMacro( NumberRows, unsigned int );

  /** Set the input image and reinitialize the list of images */
  void SetInput( const InputImageType * img );
  void SetInput( unsigned int id, const InputImageType * img );
  void AddImage( const InputImageType * img );

  const InputImageType * GetInput( void );

protected:

  ConvertImagesToCSVFilter ( void );
  ~ConvertImagesToCSVFilter( void ) {};

  virtual void GenerateData() override;

  virtual void PrintSelf( std::ostream & os, Indent indent ) const override;

private:
  ConvertImagesToCSVFilter ( const Self& );
  void operator=( const Self& );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const typename Superclass::DataObjectIdentifierType &,
    itk::DataObject * ) {};

  typename InputMaskType::Pointer                       m_InputMask;
  VnlMatrixType                                         m_VnlOutput;
  std::vector< typename InputImageType::ConstPointer >  m_ImageList;
  unsigned int                                          m_Stride;
  unsigned int                                          m_NumImages;
  unsigned int                                          m_NumberRows;

}; // End class ConvertImagesToCSVFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeConvertImagesToCSVFilter.hxx"
#endif

#endif // End !defined( _itktubeConvertImagesToCSVFilter_h )
