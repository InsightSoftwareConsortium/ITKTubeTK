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
#ifndef __itktubePadImageFilter_hxx
#define __itktubePadImageFilter_hxx

#include "itkProgressAccumulator.h"
#include "itkNumericTraits.h"
#include "itkConstantPadImageFilter.h"
#include "itkZeroFluxNeumannPadImageFilter.h"
#include "itkMirrorPadImageFilter.h"
#include "itkWrapPadImageFilter.h"
#include "itkChangeInformationImageFilter.h"

namespace itk {

namespace tube {

template <class TInputImage, class TOutputImage>
PadImageFilter<TInputImage, TOutputImage>
::PadImageFilter()
{
  m_GreatestPrimeFactor = 13;
  m_PadMethod = ZERO_FLUX_NEUMANN;
  this->SetNumberOfRequiredOutputs( 1 );
}

template <class TInputImage, class TOutputImage>
void
PadImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
 
  InputImageType * input0 = const_cast<InputImageType *>( this->GetInput(
    0 ) );
  if( !input0 )
    {
    return;
    }
 
  OutputImageType * output = this->GetOutput();
 
  RegionType region = output->GetRequestedRegion();
  region.Crop( input0->GetLargestPossibleRegion() );
  input0->SetRequestedRegion( region );
}


template <class TInputImage, class TOutputImage>
void
PadImageFilter<TInputImage, TOutputImage>
::GenerateOutputInformation()
{
  // call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();
 
  const InputImageType * input0 = this->GetInput();
  if( !input0 )
    {
    return;
    }
 
  OutputImageType * output0 = this->GetOutput();
 
  RegionType region0 = input0->GetLargestPossibleRegion();
 
  RegionType region;
  if( m_PadMethod == NO_PADDING )
    {
    region = region0;
    }
  else
    {
    // increase the size of the output by the size of the kernel
    SizeType size;
    IndexType idx;
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      long s1 = 0;
      if( m_GreatestPrimeFactor > 1 )
        {
        while( greatestPrimeFactor( region0.GetSize()[i] + s1 ) >
          m_GreatestPrimeFactor )
          {
          s1++;
          }
        }
      else if( m_GreatestPrimeFactor == 1 )
        {
        s1 += ( region0.GetSize()[i] + s1 ) % 2;
        }
      idx[i] = region0.GetIndex()[i] - s1/2;
      size[i] = region0.GetSize()[i] + s1;
      }
    region = RegionType( idx, size );
    }
  output0->SetLargestPossibleRegion( region );
}


template<class TInputImage, class TOutputImage>
void
PadImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  this->AllocateOutputs();
  const InputImageType * input0 = this->GetInput();
  OutputImageType * output0 = this->GetOutput();
  RegionType ir0 = input0->GetLargestPossibleRegion();
  RegionType or0 = output0->GetLargestPossibleRegion();

  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter( this );

  typedef typename itk::PadImageFilter
    < InputImageType, OutputImageType > PadType;
  typedef typename itk::ConstantPadImageFilter
    < InputImageType, OutputImageType > ConstantPadType;
  typedef typename itk::ZeroFluxNeumannPadImageFilter
    < InputImageType, OutputImageType > ZeroFluxPadType;
  typedef typename itk::MirrorPadImageFilter
    < InputImageType, OutputImageType > MirrorPadType;
  typedef typename itk::WrapPadImageFilter
    < InputImageType, OutputImageType > WrapPadType;
  SizeType s;
 
  typename PadType::Pointer pad0;
  switch( m_PadMethod )
    {
    case ZERO_FLUX_NEUMANN:
      {
      pad0 = ZeroFluxPadType::New();
      break;
      }
    case NO_PADDING:
    case ZERO:
      {
      pad0 = ConstantPadType::New();
      break;
      }
    case MIRROR:
      {
      pad0 = MirrorPadType::New();
      break;
      }
    case WRAP:
      {
      pad0 = WrapPadType::New();
      break;
    default:
      itkExceptionMacro( << "Unknown pad method: " << m_PadMethod );
      break;
      }
    }
  pad0->SetInput( input0 );
  pad0->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );
  if( m_PadMethod != NO_PADDING )
    {
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      s[i] = ir0.GetIndex()[i] - or0.GetIndex()[i];
      }
    pad0->SetPadLowerBound( s );
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      s[i] = or0.GetSize()[i] -
        ( ir0.GetIndex()[i] - or0.GetIndex()[i] + ir0.GetSize()[i] );
      }
    pad0->SetPadUpperBound( s );
    }
  progress->RegisterInternalFilter( pad0, 0.5f );
  pad0->GraftOutput( output0 );
  pad0->Update();
  this->GraftOutput( pad0->GetOutput() );
}


template<class TInputImage, class TOutputImage>
void
PadImageFilter<TInputImage, TOutputImage>
::PrintSelf( std::ostream &os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "GreatestPrimeFactor: "  << m_GreatestPrimeFactor
    << std::endl;
  os << indent << "PadMethod: "  << m_PadMethod << std::endl;
}
 
}// end namespace tube

}// end namespace itk
#endif
