/*=========================================================================

Library:   TubeTK/VTree

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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
#ifndef __itkTubeExtractor_txx
#define __itkTubeExtractor_txx

#include "itkTubeExtractor.h"

namespace itk
{
 
/**
 * Constructor */
template<class TInputImage>
TubeExtractor<TInputImage>
::TubeExtractor()
{
  m_Debug = true;
  m_IdleCallBack = NULL;
  m_StatusCallBack = NULL;
  m_NewTubeCallBack = NULL;
  m_AbortProcess = NULL;
}

/**
 * Destructor */
template<class TInputImage>
TubeExtractor<TInputImage>
::~TubeExtractor()
{
  //delete m_X;
}

/**
 * Set the input image */
template<class TInputImage>
void
TubeExtractor<TInputImage>
::SetInputImage(ImagePointer inputImage )
{
  m_Image = inputImage;

  m_RidgeOp = RidgeExtractor<ImageType>::New();
  m_RidgeOp->SetInputImage(m_Image);
 
  m_RadiusOp = RadiusExtractor<ImageType>::New();
  m_RadiusOp->SetInputImage(m_Image);
  m_RidgeOp->SetRadiusExtractor(m_RadiusOp);
}

/**
 * Set Data Min value */
template<class TInputImage>
void
TubeExtractor<TInputImage>
::SetDataMin(double dataMin)
{
  m_RidgeOp->SetDataMin(dataMin);
  m_RadiusOp->SetDataMin(dataMin); 
}

/**
 * Get Data Min value */
template<class TInputImage>
double
TubeExtractor<TInputImage>
::GetDataMin(void)
{
  return m_RidgeOp->GetDataMin(); 
}

/**
 * Set Data Max value */
template<class TInputImage>
void
TubeExtractor<TInputImage>
::SetDataMax(double dataMax)
{
  m_RidgeOp->SetDataMin(dataMax);
  m_RadiusOp->SetDataMin(dataMax); 
}


/**
 * Get Data Max value */
template<class TInputImage>
double
TubeExtractor<TInputImage>
::GetDataMax(void)
{
  return m_RidgeOp->GetDataMax(); 
}

/**
 * Set Radius*/
template<class TInputImage>
void
TubeExtractor<TInputImage>
::SetRadius(double radius)
{
  m_RidgeOp->SetScale(radius);
  m_RadiusOp->SetRadius0(radius);
}

/**
 * Get Radius*/
template<class TInputImage>
double
TubeExtractor<TInputImage>
::GetRadius(void)
{
  return m_RidgeOp->GetScale();
}

/**
 * Extract the ridge */
template<class TInputImage>
void
TubeExtractor<TInputImage>
::ExtractRidge(bool extractRidge)
{
  m_RadiusOp->SetExtractRidge(extractRidge);
}

/**
 * Get Extract ridge */
template<class TInputImage>
bool
TubeExtractor<TInputImage>
::ExtractRidge(void)
{
  return m_RidgeOp->GetExtractRidge();
}

/**
 * Extract a valley */
template<class TInputImage>
void
TubeExtractor<TInputImage>
::ExtractValley(bool extractValley)
{
  //m_RidgeOp->SetExtractValley(extractValley);
  m_RadiusOp->SetExtractRidge(!extractValley);
}

/**
 * Get Extract valley */
template<class TInputImage>
bool
TubeExtractor<TInputImage>
::ExtractValley(void)
{
  return m_RidgeOp->GetExtractValley();
}

/**
 * Get the ridge extractor */
template<class TInputImage>
typename RidgeExtractor<TInputImage>::Pointer
TubeExtractor<TInputImage>
::GetRidgeOp(void)
{
  return &  m_RidgeOp;
}
  
/**
 * Get the radius extractor */
template<class TInputImage>
typename RadiusExtractor<TInputImage>::Pointer
TubeExtractor<TInputImage>
::GetRadiusOp(void)
{
  return & m_RadiusOp;
}
  
/**
 * Extract the tube given the position of the first point
 * and the tube ID */
template<class TInputImage> 
bool 
TubeExtractor<TInputImage>
::ExtractTube(float x, float y, float z, unsigned int tubeID)
{

  m_X[0] = x;
  m_X[1] = y;
  m_X[2] = z;

  m_Tube = m_RidgeOp->Extract(m_X, tubeID);
  
  if(!m_Tube)
    {
    if(m_Debug) std::cout << "m_RidgeOp->Extract() fails !" << std::endl;
    return false;
    }

  // Set the Spacing of the tube as the same spacing of the image
  double spacing[3];
  for(unsigned int i=0;i<3;i++)
    {
    spacing[i] = m_Image->GetSpacing()[i];
    }
  m_Tube->GetIndexToObjectTransform()->SetScaleComponent(spacing);

  m_RidgeOp->SmoothTubeX(m_Tube, 15);

  if(m_AbortProcess != NULL)
    {
    if(m_AbortProcess()) 
      {
      if(m_StatusCallBack) 
        { 
        m_StatusCallBack("Extract: Ridge", "Aborted", 0);
        }
      return false;
      }
    } 

  double tR = m_RadiusOp->GetRadius0();

  if(m_RidgeOp->GetDynamicScale())
    {
    m_RadiusOp->SetRadius0(m_RidgeOp->GetDynamicScaleUsed());
    }

  if(!m_RadiusOp->CalcRadii(m_Tube)) 
    {
    m_RadiusOp->SetRadius0(tR);
    return false;
    }
  m_RadiusOp->SetRadius0(tR);

  if(m_NewTubeCallBack != NULL)
    {
    m_NewTubeCallBack(m_Tube);
    }

  if(m_StatusCallBack) 
    {
    char s[80];
    sprintf(s, "%d points", m_Tube->GetPoints().size());
    m_StatusCallBack("Extract: Ridge", s, 0);
    }

  return true;

}

/**
 * Delete a tube */
template<class TInputImage> 
bool 
  TubeExtractor<TInputImage>
::DeleteTube(TubeType* tube)
{
  return this->M_RidgeOp->deleteTube(tube);
}

/** 
 * Get the last tube extracted */
template<class TInputImage> 
typename TubeExtractor<TInputImage>::TubePointer
  TubeExtractor<TInputImage>
::GetLastTube(void)
{
  return m_Tube;
}

/**
 * Set the tube color */
template<class TInputImage> 
void
  TubeExtractor<TInputImage>
::SetColor(float* color)
{
  m_Color = color;
}


/**
 * Get the last position */
template<class TInputImage> 
typename TubeExtractor<TInputImage>::ContinuousIndexType
  TubeExtractor<TInputImage>
::GetLastPosition(void)
{
  return m_X;
}

/**
 * Set the idle call back */
template<class TInputImage>
void
  TubeExtractor<TInputImage>
::IdleCallBack(bool (*idleCallBack)())
{
  m_IdleCallBack = idleCallBack;
}

/**
 * Set the status callback  */
template<class TInputImage>
void
  TubeExtractor<TInputImage>
::StatusCallBack(void (*statusCallBack)(char *, char *, int))
{
  m_StatusCallBack = statusCallBack;
  m_RidgeOp->StatusCallBack(statusCallBack);
  m_RadiusOp->StatusCallBack(statusCallBack);
}

/**
 * Set the status callback  */
template<class TInputImage>
void
  TubeExtractor<TInputImage>
::NewTubeCallBack(void (*newTubeCallBack)(TubeType* ))
{
  m_NewTubeCallBack = newTubeCallBack;
}

/**
 * Abort the process  */
template<class TInputImage>
void
  TubeExtractor<TInputImage>
::AbortProcess(bool (*abortProcess)())
{
  m_AbortProcess = abortProcess;
}


}; // end namespace itk

#endif
