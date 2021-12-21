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
#ifndef __tubeEnhanceTubesUsingDiffusion_hxx
#define __tubeEnhanceTubesUsingDiffusion_hxx


#include <vector>

namespace tube
{

template< class TPixel, unsigned int Dimension >
EnhanceTubesUsingDiffusion< TPixel, Dimension >
::EnhanceTubesUsingDiffusion( void )
{
  m_MinSigma = 1.0;
  m_MaxSigma = 1.0;
  m_NumSigmaSteps = 1;

  m_Filter = FilterType::New();
  m_Filter->SetDefaultPars();
}

template< class TPixel, unsigned int Dimension >
void
EnhanceTubesUsingDiffusion< TPixel, Dimension >
::Update( void )
{
  // set scales
  std::vector< float > scales( m_NumSigmaSteps );

  double deltaSigma = m_MaxSigma - m_MinSigma;

  for( unsigned int i = 0; i < m_NumSigmaSteps; i++ )
    {
    scales[i] = m_MinSigma + i * ( deltaSigma / m_NumSigmaSteps );
    }

  m_Filter->SetScales( scales );

  // compute vesselness image
  m_Filter->Update();
}

template< class TPixel, unsigned int Dimension >
void
EnhanceTubesUsingDiffusion< TPixel, Dimension >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "MinSigma             : " << GetMinSigma() << std::endl;
  os << indent << "MaxSigma             : " << GetMaxSigma() << std::endl;
  os << indent << "NumSigmaSteps        : " << GetNumSigmaSteps() << std::endl;

  os << indent << "RecalculateTubeness  : " << GetRecalculateTubeness()
                                            << std::endl;
  os << indent << "Beta                 : " << GetBeta() << std::endl;
  os << indent << "Gamma                : " << GetGamma() << std::endl;
  os << indent << "Epsilon              : " << GetEpsilon() << std::endl;
  os << indent << "Omega                : " << GetOmega() << std::endl;
  os << indent << "Sensitivity          : " << GetSensitivity() << std::endl;

  os << indent << "TimeStep             : " << GetTimeStep()  << std::endl;
  os << indent << "Iterations           : " << GetIterations() << std::endl;
}

}

#endif
