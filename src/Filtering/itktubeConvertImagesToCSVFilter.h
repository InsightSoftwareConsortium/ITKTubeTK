/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itktubeConvertImagesToCSVFilter_h
#define itktubeConvertImagesToCSVFilter_h

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

template <class TInputImage, class TInputMask>
class ConvertImagesToCSVFilter : public ProcessObject
{
public:
  /** Standard class type alias. */
  using Self = ConvertImagesToCSVFilter;
  using Superclass = ProcessObject;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using InputImageType = TInputImage;
  using ImagePointer = typename InputImageType::Pointer;
  using InputPixelType = typename InputImageType::PixelType;
  using IndexType = typename InputImageType::IndexType;
  using ReaderType = ImageFileReader<InputImageType>;
  using InputImageIteratorType = ImageRegionConstIterator<InputImageType>;

  using InputMaskType = TInputMask;
  using MaskPointer = typename InputMaskType::Pointer;
  using MaskPixelType = typename InputMaskType::PixelType;
  using MaskIndexType = typename InputMaskType::IndexType;
  using MaskReaderType = ImageFileReader<InputMaskType>;
  using MaskIteratorType = ImageRegionConstIterator<InputMaskType>;

  using VnlMatrixType = vnl_matrix<InputPixelType>;
  using OutputType = SimpleDataObjectDecorator<VnlMatrixType>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro(ConvertImagesToCSVFilter);

  /** ImageDimension constants */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputMask::ImageDimension);

  itkSetObjectMacro(InputMask, InputMaskType);
  itkGetModifiableObjectMacro(InputMask, InputMaskType);

  OutputType *
  GetOutput();

  itkGetMacro(Stride, unsigned int);

  itkSetClampMacro(Stride, unsigned int, 1, std::numeric_limits<unsigned int>::max());

  itkGetMacro(NumImages, unsigned int);
  itkSetMacro(NumImages, unsigned int);

  itkSetMacro(NumberRows, unsigned int);
  itkGetMacro(NumberRows, unsigned int);

  /** Set the input image and reinitialize the list of images */
  void
  SetInput(const InputImageType * img);
  void
  SetInput(unsigned int id, const InputImageType * img);
  void
  AddImage(const InputImageType * img);

  const InputImageType *
  GetInput(void);

protected:
  ConvertImagesToCSVFilter(void);
  ~ConvertImagesToCSVFilter(void) {};

  virtual void
  GenerateData() override;

  virtual void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  ConvertImagesToCSVFilter(const Self &);
  void
  operator=(const Self &);

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void
  SetInput(const typename Superclass::DataObjectIdentifierType &, itk::DataObject *) override {};

  typename InputMaskType::Pointer                    m_InputMask;
  VnlMatrixType                                      m_VnlOutput;
  std::vector<typename InputImageType::ConstPointer> m_ImageList;
  unsigned int                                       m_Stride;
  unsigned int                                       m_NumImages;
  unsigned int                                       m_NumberRows;

}; // End class ConvertImagesToCSVFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeConvertImagesToCSVFilter.hxx"
#endif

#endif // End !defined( _itktubeConvertImagesToCSVFilter_h )
