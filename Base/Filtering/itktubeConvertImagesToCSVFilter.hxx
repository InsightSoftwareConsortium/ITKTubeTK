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

#ifndef __itktubeConvertImagesToCSVFilter_hxx
#define __itktubeConvertImagesToCSVFilter_hxx

#include "itktubeConvertImagesToCSVFilter.h"

namespace itk
{

namespace tube
{

/** Constructor */
template< class TInputImage, class TInputMask >
ConvertImagesToCSVFilter< TInputImage, TInputMask >
::ConvertImagesToCSVFilter( void )
{
  m_InputMask = NULL;
  m_NumberRows = 0;
  this->ProcessObject::SetNthOutput( 0, OutputType::New().GetPointer() );
}

template< class TInputImage, class TInputMask >
void
ConvertImagesToCSVFilter< TInputImage, TInputMask >
::GenerateData(void)
{
  const unsigned int ARows =
    m_InputMask->GetLargestPossibleRegion().GetNumberOfPixels() / m_Stride;
  const unsigned int ACols = m_ImageList.size() + 1;
  m_VnlOutput.set_size( ARows, ACols );
  std::vector< IteratorType * > iterList;
  for( unsigned int i = 0; i < m_NumImages; ++i )
    {
    iterList.push_back( new IteratorType(m_ImageList[i],
    m_ImageList[i]->GetLargestPossibleRegion()) );
    }
  MaskIteratorType maskIter( m_InputMask,
    m_InputMask->GetLargestPossibleRegion() );
  unsigned int i = 0;
  while( !maskIter.IsAtEnd() )
    {
    if( maskIter.Get() != 0 )
      {
      for( i = 0; i<m_NumImages; ++i )
        {
        m_VnlOutput( m_NumberRows, i ) = iterList[i]->Get();
        }
      m_VnlOutput( m_NumberRows, i ) = maskIter.Get();
      m_NumberRows++;
      }
    for( int s = 0; s<m_Stride && !maskIter.IsAtEnd(); ++s )
      {
      for( unsigned int i = 0; i<m_NumImages; ++i )
        {
        ++( *iterList[i] );
        }
      ++maskIter;
      }
    }
  for( unsigned int i = 0; i<iterList.size(); ++i )
    {
    delete iterList[i];
    }
  iterList.clear();
  typename OutputType::Pointer outputPtr = this->GetOutput();
  outputPtr->Set( m_VnlOutput );
}

template< class TInputImage, class TInputMask >
SimpleDataObjectDecorator<vnl_matrix <typename TInputImage::PixelType> >*
ConvertImagesToCSVFilter< TInputImage, TInputMask >
::GetOutput()
{
  // we assume that the first output is of the templated type
  return itkDynamicCastInDebugMode< OutputType * >( this->GetPrimaryOutput() );
}

/** Set the input image and reinitialize the list of images */
template< class TInputImage, class TInputMask >
void
ConvertImagesToCSVFilter< TInputImage, TInputMask >
::SetInput(InputImageType* image)
{
  m_ImageList.clear();
  this->m_ImageList.push_back( image );
  this->Modified();
}

template< class TInputImage, class TInputMask >
const TInputImage*
ConvertImagesToCSVFilter< TInputImage, TInputMask >
::GetInput()
{
  return *m_ImageList.begin();
}

template< class TInputImage, class TInputMask >
void
ConvertImagesToCSVFilter< TInputImage, TInputMask >
::AddImage(InputImageType* image)
{
  this->m_ImageList.push_back( image );
  this->Modified();
}

/** PrintSelf */
template< class TInputImage, class TInputMask >
void
ConvertImagesToCSVFilter< TInputImage, TInputMask >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "ImageList size = " << m_ImageList.size() << std::endl;
  if (m_ImageList.size() > 0)
  {
    os << indent << "ImageList[0] = " << m_ImageList[0] << std::endl;
  }
  os << indent << "InputMask = " << m_InputMask << std::endl;
  os << indent << "VnlOutput = " << m_VnlOutput << std::endl;
  os << indent << "Stride = " << m_Stride << std::endl;
  os << indent << "NumImages = " << m_NumImages << std::endl;
  os << indent << "NumberRows = " << m_NumberRows << std::endl;
}

} // End namespace tube

} // End namespace itk
#endif // End !defined( __itktubeConvertImagesToCSVFilter_hxx )
