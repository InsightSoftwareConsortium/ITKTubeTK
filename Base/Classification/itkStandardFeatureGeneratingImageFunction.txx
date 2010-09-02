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
#ifndef __itkStandardFeatureGeneratingImageFunction_txx
#define __itkStandardFeatureGeneratingImageFunction_txx

#include "itkStandardFeatureGeneratingImageFunction.h"

#include "itkRidgeExtractor.h"
#include "itkJointHistogramImageFunction.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkDivideImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNeighborhoodIterator.h"
#include "vnl/vnl_math.h"
#include "math.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>

namespace itk
{

/**
 * Set the input Image
 */
template <class TInputImage, class TCoordRep>
StandardFeatureGeneratingImageFunction<TInputImage,TCoordRep>
::StandardFeatureGeneratingImageFunction()
{
  this->m_Features.push_back("RidgenessDiffSmall");
  this->m_Features.push_back("RidgenessDiffMedium");
  this->m_Features.push_back("RidgenessDiffLarge");

  this->m_Features.push_back("NormalizedScaledIntensityDiff");
  this->m_Features.push_back("ScaledIntensityDiff");

  this->m_Features.push_back("ZScoreAdds");
  this->m_Features.push_back("ZScoreSubs");
  this->m_Features.push_back("ZScoreNoms");
}

template <class TInputImage, class TCoordRep>
void
StandardFeatureGeneratingImageFunction<TInputImage,TCoordRep>
::PrepFilter()
{
  m_InputCalcSmall = CalculatorType::New();
  m_InputCalcSmall->SetInputImage( this->GetInputImage() );
  m_InputDataMin = m_InputCalcSmall->GetDataMin();
  m_InputDataMax = m_InputCalcSmall->GetDataMin();
  m_InputCalcSmall->SetDataMin( m_InputDataMax );
  m_InputCalcSmall->SetDataMax( m_InputDataMin );

  m_InputCalcMedium = CalculatorType::New();
  m_InputCalcMedium->SetInputImage( this->GetInputImage() );
  m_InputCalcMedium->SetDataMin( m_InputDataMax );
  m_InputCalcMedium->SetDataMax( m_InputDataMin );

  m_InputCalcLarge = CalculatorType::New();
  m_InputCalcLarge->SetInputImage( this->GetInputImage() );
  m_InputCalcLarge->SetDataMin( m_InputDataMax );
  m_InputCalcLarge->SetDataMax( m_InputDataMin );

  m_PriorCalcSmall = CalculatorType::New();
  m_PriorCalcSmall->SetInputImage( m_Prior );
  m_PriorDataMax = m_PriorCalcSmall->GetDataMax();
  m_PriorDataMin = m_PriorCalcSmall->GetDataMin();
  m_PriorCalcSmall->SetDataMin( m_PriorDataMax );
  m_PriorCalcSmall->SetDataMax( m_PriorDataMin );

  m_PriorCalcMedium = CalculatorType::New();
  m_PriorCalcMedium->SetInputImage( m_Prior );
  m_PriorCalcMedium->SetDataMin( m_PriorDataMax );
  m_PriorCalcMedium->SetDataMax( m_PriorDataMin );

  m_PriorCalcLarge = CalculatorType::New();
  m_PriorCalcLarge->SetInputImage( m_Prior );
  m_PriorCalcLarge->SetDataMin( m_PriorDataMax );
  m_PriorCalcLarge->SetDataMax( m_PriorDataMin );

  m_AddJHCalc = HistCalcType::New();
  m_AddJHCalc->SetInputImage( this->GetInputImage() ); 
  m_AddJHCalc->SetInputMask( m_Prior );
  m_AddJHCalc->SetMeanHistogram( m_AddMeanHist );
  m_AddJHCalc->SetStandardDeviationHistogram( m_AddStdevHist );

  m_SubJHCalc = HistCalcType::New();
  m_SubJHCalc->SetInputImage( this->GetInputImage() ); 
  m_SubJHCalc->SetInputMask( m_Prior );
  m_SubJHCalc->SetMeanHistogram( m_SubMeanHist );
  m_SubJHCalc->SetStandardDeviationHistogram( m_SubStdevHist );

  m_NomJHCalc = HistCalcType::New();
  m_NomJHCalc->SetInputImage( this->GetInputImage() ); 
  m_NomJHCalc->SetInputMask( m_Prior );
  m_NomJHCalc->SetMeanHistogram( m_NormMeanHist );
  m_NomJHCalc->SetStandardDeviationHistogram( m_NormStdevHist );

  m_InputCalcSmall->SetScale( m_SigmaSmall );
  m_InputCalcMedium->SetScale( m_SigmaMedium );
  m_InputCalcLarge->SetScale( m_SigmaLarge );

  m_PriorCalcSmall->SetScale( m_SigmaSmall );
  m_PriorCalcMedium->SetScale( m_SigmaMedium );
  m_PriorCalcLarge->SetScale( m_SigmaLarge );
}

template <class TInputImage, class TCoordRep>
typename StandardFeatureGeneratingImageFunction<TInputImage,TCoordRep>::OutputType
StandardFeatureGeneratingImageFunction<TInputImage,TCoordRep>
::EvaluateAtIndex( const IndexType & index ) const
{
  PointType curPoint;
  this->GetInputImage()->TransformIndexToPhysicalPoint( index, curPoint );

  // holders for features
  double v1sg, v1mg, v1lg, v2g, v3g, v1se, v1me, v1le, v2e, v3e;
  double roundness = 0;
  double curvature = 0;
  double zAdd, zSub, zNom;
      
  v1sg = m_PriorCalcSmall->Ridgeness( index, roundness,
                                    curvature );
  v1mg = m_PriorCalcMedium->Ridgeness( index, roundness,
                                     curvature );
  v1lg = m_PriorCalcLarge->Ridgeness( index, roundness,
                                    curvature );

  v1se = m_InputCalcSmall->Ridgeness( index, roundness,
                                    curvature );
  v1me = m_InputCalcMedium->Ridgeness( index, roundness,
                                     curvature );
  v1le = m_InputCalcLarge->Ridgeness( index, roundness,
                                    curvature );
  
  double pS = m_PriorCalcSmall->Intensity( index );
  double pM = m_PriorCalcMedium->Intensity( index );
  double pL = m_PriorCalcLarge->Intensity( index );

  double iS = m_InputCalcSmall->Intensity( index );
  double iM = m_InputCalcMedium->Intensity( index );
  double iL = m_InputCalcLarge->Intensity( index );

  if( pL != 0 )
    {
    v2g = ( pS - pL ) / pL;
    }
  else
    {
    v2g = 0;
    }
  if( iL != 0 )
    {
    v2e = ( iS - iL ) / iL;
    }
  else
    {
    v2e = 0;
    }
      
  v3g = pM;
  v3e = iM;
      
  zAdd = m_AddJHCalc->Evaluate( curPoint );
  zSub = m_SubJHCalc->Evaluate( curPoint );
  zNom = m_NomJHCalc->Evaluate( curPoint );

  OutputType out( this->m_Features.size() );

  out[0] = curPoint[0];
  out[1] = curPoint[1];
  out[2] = v1sg-v1se;
  out[3] = v1mg-v1me;
  out[4] = v1lg-v1le; 
  out[5] = v2g-v2e;
  out[6] = v3g-v3e;
  out[7] = zAdd;
  out[8] = zSub;
  out[9] = zNom;

  return out;
}

}

#endif
