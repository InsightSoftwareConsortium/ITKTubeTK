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
#ifndef __itktubeComputeTrainingMask_hxx
#define __itktubeComputeTrainingMask_hxx

#include "itktubeComputeTrainingMask.h"

namespace itk
{
namespace tube
{

template< class TInputImage >
ComputeTrainingMask< TInputImage >
::ComputeTrainingMask()
{
  m_Gap = 0;
  m_NotVesselWidth = 1.0;
  m_BinaryThinning = BinaryThinningFilterType::New();

  m_Threshold = ThresholdFilterType::New();

  m_Threshold->SetLowerThreshold( 0);
  m_Threshold->SetInsideValue( 0 );
  m_Threshold->SetOutsideValue( 255 );

  BallType m_Ball;
  m_Ball.SetRadius( 1 );
  m_Ball.CreateStructuringElement();

  m_Dilate = DilateFilterType::New();
  m_Dilate->SetObjectValue(255 );
  m_Dilate->SetKernel( m_Ball );

  m_Substract = SubstractFilterType::New();
  m_MultiplyCenterLine = MultiplyFilterType::New();
  m_MultiplyCenterLine->SetConstant(255);

  m_DivideImage = DivideFilterType::New();
  m_DivideImage->SetConstant(2.0);

  m_Add = AddFilterType::New();
  m_Cast = CastFilterType::New();
  m_CastNotVessel = CastFilterType::New();

  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfRequiredOutputs(2);

  typename ImageTypeShort::Pointer output1 =
    static_cast< ImageTypeShort * >( this->MakeOutput(1).GetPointer() );
  this->ProcessObject::SetNthOutput( 1,output1);
}

template< class TInputImage >
void
ComputeTrainingMask< TInputImage >
::ApplyDilateMorphologyFilter( typename ImageType::Pointer &input )
{
  for ( int r = 0; r<m_NotVesselWidth; r++ )
    {
    m_Dilate->SetInput(input);
    m_Dilate->Update();
    input=m_Dilate->GetOutput();
    input->DisconnectPipeline();
    }
  return;
}

template< class TInputImage >
void
ComputeTrainingMask< TInputImage >
::GenerateData()
{
  typename ImageType::Pointer input = ImageType::New();
  input->Graft( const_cast< ImageType * >( this->GetInput() ));
  m_BinaryThinning->SetInput(input);
  m_Threshold->SetInput( input );
  m_Threshold->SetUpperThreshold( m_Gap );
  m_Threshold->Update();
  typename ImageType::Pointer image = m_Threshold->GetOutput();

  ApplyDilateMorphologyFilter( image );
  typename ImageType::Pointer dilatedImage = image;
  ApplyDilateMorphologyFilter( image );
  m_Substract->SetInput1(image);
  m_Substract->SetInput2(dilatedImage);

  m_MultiplyCenterLine->SetInput(m_BinaryThinning->GetOutput());
  m_DivideImage->SetInput(m_Substract->GetOutput());
  m_Add->SetInput1(m_MultiplyCenterLine->GetOutput());
  m_Add->SetInput2(m_DivideImage->GetOutput());

  m_CastNotVessel->SetInput(m_Substract->GetOutput());
  m_CastNotVessel->GraftOutput(const_cast< ImageTypeShort * >(this->GetOutput(1)));
  m_CastNotVessel->Update();
  this->GraftNthOutput(1,m_CastNotVessel->GetOutput());

  m_Cast->SetInput(m_Add->GetOutput());
  m_Cast->GraftOutput(this->GetOutput());
  m_Cast->Update();
  this->GraftOutput( m_Cast->GetOutput() );
}

template< class TInputImage >
const typename ComputeTrainingMask< TInputImage >::ImageTypeShort*
ComputeTrainingMask< TInputImage >
::GetNotVesselMask()
{
  return itkDynamicCastInDebugMode< ImageTypeShort * >( this->GetOutput(1) );
}

template< class TInputImage >
ComputeTrainingMask< TInputImage >
::~ComputeTrainingMask()
{
}

template< class TInputImage >
void
ComputeTrainingMask< TInputImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  os << indent << "Gap = " << m_Gap << std::endl;
  os << indent << "NotVesselWidth = " << m_NotVesselWidth << std::endl;
  ImageType* inputPtr = const_cast< ImageType * >( this->GetInput() );
  if( inputPtr )
    {
    os << indent << "Input = " << inputPtr << std::endl;
    }
}

}//end of tube namespace
}// end of itk namespace
#endif
