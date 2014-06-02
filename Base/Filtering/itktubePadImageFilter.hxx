/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itktubePadImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itktubePadImageFilter_txx
#define __itktubePadImageFilter_txx

#include "itktubePadImageFilter.h"
#include "itkProgressAccumulator.h"
#include "itkNumericTraits.h"
#include "itkConstantPadImageFilter.h"
#include "itkZeroFluxNeumannPadImageFilter.h"
#include "itkMirrorPadImageFilter.h"
#include "itkWrapPadImageFilter.h"
#include "itkChangeInformationImageFilter.h"

namespace itk {

namespace tube {

template <class TInputImage, class TInputKernel, class TOutputImage, class TKernelOutput>
PadImageFilter<TInputImage, TInputKernel, TOutputImage, TKernelOutput>
::PadImageFilter()
{
  m_GreatestPrimeFactor = 13;
  m_PadMethod = ZERO_FLUX_NEUMANN;
//  this->SetNumberOfRequiredInputs(2);
  this->SetNumberOfRequiredOutputs(2);
  this->SetNthOutput( 1, TKernelOutput::New() );
}

template <class TInputImage, class TInputKernel, class TOutputImage, class TKernelOutput>
void 
PadImageFilter<TInputImage, TInputKernel, TOutputImage, TKernelOutput>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  InputImageType * input0 = const_cast<InputImageType *>(this->GetInput(0));
  InputImageType * input1 = const_cast<InputImageType *>(this->GetInput(1));
  if ( !input0 )
    { 
    return;
    }
  
  OutputImageType * output = this->GetOutput();
  
  RegionType region = output->GetRequestedRegion();
  region.Crop( input0->GetLargestPossibleRegion() );
  input0->SetRequestedRegion( region );
  
  // input1 is not required
  if ( !input1 )
    { 
    return;
    }
  
  region = output->GetRequestedRegion();
  region.Crop( input1->GetLargestPossibleRegion() );
  input1->SetRequestedRegion( region );
}


template <class TInputImage, class TInputKernel, class TOutputImage, class TKernelOutput>
void 
PadImageFilter<TInputImage, TInputKernel, TOutputImage, TKernelOutput>
::GenerateOutputInformation()
{
  // call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();
  
  const InputImageType * input0 = this->GetInput();
  const InputKernelType * input1 = this->GetInputKernel();
  if ( !input0 )
    { 
    return;
    }
  
  OutputImageType * output0 = this->GetOutput();
  OutputKernelType * output1 = this->GetOutputKernel();
  
  RegionType region0 = input0->GetLargestPossibleRegion();
  RegionType region1; // set later
  
  // dummy region to avoid code duplication when there is no input1
  IndexType nullidx;
  nullidx.Fill(0);
  SizeType nullsize;
  nullsize.Fill(0);
  RegionType nullregion( nullidx, nullsize );
    
  if( input1 )
    {
    region1 = input1->GetLargestPossibleRegion();
    }
  else
    {
    region1 = nullregion;
    }
  
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
    for( int i=0; i<ImageDimension; i++ )
      {
      long s1 = std::max( (long)region1.GetSize()[i] - 1, (long)0 );
      if( m_GreatestPrimeFactor > 1 )
        {
        while( greatestPrimeFactor( region0.GetSize()[i] + s1 ) > m_GreatestPrimeFactor )
          {
          s1++;
          }
        }
      else if( m_GreatestPrimeFactor == 1 )
        {
        // make sure the total size is even
        s1 += ( region0.GetSize()[i] + s1 ) % 2;
        }
      idx[i] = region0.GetIndex()[i] - s1/2;
      size[i] = region0.GetSize()[i] + s1;
      }
    region = RegionType( idx, size );
    }
  output0->SetLargestPossibleRegion( region );
  // make sure that output1 is actually there - it can be set to NULL by subclasses
  if( output1 )
    {
    if( input1 )
      {
      output1->SetLargestPossibleRegion( region );
      }
    else
      {
      output1->SetLargestPossibleRegion( nullregion );
      }
    }
  // std::cout << region << std::endl;
}


template<class TInputImage, class TInputKernel, class TOutputImage, class TKernelOutput>
void
PadImageFilter<TInputImage, TInputKernel, TOutputImage, TKernelOutput>
::GenerateData()
{
  this->AllocateOutputs();
  const InputImageType * input0 = this->GetInput();
  const InputKernelType * input1 = this->GetInputKernel();
  OutputImageType * output0 = this->GetOutput();
  OutputKernelType * output1 = this->GetOutputKernel();
  RegionType ir0 = input0->GetLargestPossibleRegion();
  RegionType ir1;  // set later
  RegionType or0 = output0->GetLargestPossibleRegion();
  RegionType or1 = output1->GetLargestPossibleRegion();

  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  typedef typename itk::PadImageFilter< InputImageType, OutputImageType > PadType;
  typedef typename itk::ConstantPadImageFilter< InputImageType, OutputImageType > ConstantPadType;
  typedef typename itk::ZeroFluxNeumannPadImageFilter< InputImageType, OutputImageType > ZeroFluxPadType;
  typedef typename itk::MirrorPadImageFilter< InputImageType, OutputImageType > MirrorPadType;
  typedef typename itk::WrapPadImageFilter< InputImageType, OutputImageType > WrapPadType;
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
      itkExceptionMacro(<< "Unknown pad method: " << m_PadMethod);
      break;
      }
    }
  pad0->SetInput( input0 );
  pad0->SetNumberOfThreads( this->GetNumberOfThreads() );
  if( m_PadMethod != NO_PADDING )
    {
    for( int i=0; i<ImageDimension; i++ )
      {
      s[i] = ir0.GetIndex()[i] - or0.GetIndex()[i];
      }
    pad0->SetPadLowerBound( s );
    for( int i=0; i<ImageDimension; i++ )
      {
      s[i] = or0.GetSize()[i] - ( ir0.GetIndex()[i] - or0.GetIndex()[i] + ir0.GetSize()[i]);
      }
    pad0->SetPadUpperBound( s );
    }
  progress->RegisterInternalFilter( pad0, 0.5f );
  pad0->GraftOutput( output0 );
  pad0->Update();
  this->GraftOutput( pad0->GetOutput() );


  if( input1 )
    {
    ir1 = input1->GetLargestPossibleRegion();
  
    typedef typename itk::ConstantPadImageFilter< InputKernelType, OutputKernelType > KernelPadType;
    typename KernelPadType::Pointer pad1 = KernelPadType::New();
    pad1->SetInput( input1 );
    pad1->SetNumberOfThreads( this->GetNumberOfThreads() );
    for( int i=0; i<ImageDimension; i++ )
      {
      s[i] = ( or1.GetSize()[i] - ir1.GetSize()[i] ) / 2;
      }
    pad1->SetPadUpperBound( s );
    for( int i=0; i<ImageDimension; i++ )
      {
      // s[i] = itk::Math::Ceil(( or1.GetSize()[i] - ir1.GetSize()[i] ) / 2.0 );
      // this line should do the same, but without requirement on ITK cvs
      s[i] = ( or1.GetSize()[i] - ir1.GetSize()[i] ) / 2 +  ( or1.GetSize()[i] - ir1.GetSize()[i] ) % 2;
      }
    pad1->SetPadLowerBound( s );
    progress->RegisterInternalFilter( pad1, 0.5f );
    
    typedef typename itk::ChangeInformationImageFilter< OutputKernelType > ChangeType;
    typename ChangeType::Pointer change = ChangeType::New();
    change->SetInput( pad1->GetOutput() );
    change->SetUseReferenceImage( true );
    change->SetReferenceImage( output1 );
    change->SetChangeRegion( true );
    // no progress for change - it does almost nothing
    
    change->GraftOutput( output1 );
    change->Update();
    this->GraftNthOutput( 1, change->GetOutput() );
    }
}


template<class TInputImage, class TInputKernel, class TOutputImage, class TKernelOutput>
void
PadImageFilter<TInputImage, TInputKernel, TOutputImage, TKernelOutput>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "GreatestPrimeFactor: "  << m_GreatestPrimeFactor << std::endl;
  os << indent << "PadMethod: "  << m_PadMethod << std::endl;
}
  
}// end namespace tube

}// end namespace itk
#endif
