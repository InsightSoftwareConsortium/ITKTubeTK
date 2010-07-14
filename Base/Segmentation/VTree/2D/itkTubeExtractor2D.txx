/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkTubeExtractor2D.txx,v $
  Language:  C++
  Date:      $Date: 2005/10/18 00:52:53 $
  Version:   $Revision: 1.19 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTubeExtractor2D_txx
#define __itkTubeExtractor2D_txx

#include "itkTubeExtractor2D.h"

namespace itk
{
 
/**
 * Constructor */
template<class TInputImage>
TubeExtractor2D<TInputImage>
::TubeExtractor2D()
{
  m_Debug = true;
  m_IdleCallBack = NULL;
  m_StatusCallBack = NULL;
  m_NewTubeCallBack = NULL;
  m_AbortProcess = NULL;
  m_X = new VectorType(ImageDimension);
}

/**
 * Destructor */
template<class TInputImage>
TubeExtractor2D<TInputImage>
::~TubeExtractor2D()
{
  delete m_X;
}

/**
 * Set the input image */
template<class TInputImage>
void
TubeExtractor2D<TInputImage>
::SetInputImage(ImagePointer inputImage )
{

  m_Image = inputImage;
  this->m_RidgeOp = RidgeExtractor2D<ImageType>::New();
  
  this->m_RidgeOp->SetInputImage(m_Image);
  this->m_RadiusOp = RadiusExtractor2D<ImageType>::New();
  this->m_RadiusOp->SetInputImage(m_Image);
  this->m_RidgeOp->SetRadiusExtractor(this->m_RadiusOp);
}

/**
 * Set Data Min value */
template<class TInputImage>
void
TubeExtractor2D<TInputImage>
::SetDataMin(double dataMin)
{
  this->m_RidgeOp->SetDataMin(dataMin);
  this->m_RadiusOp->SetDataMin(dataMin); 
}

/**
 * Get Data Min value */
template<class TInputImage>
double
TubeExtractor2D<TInputImage>
::GetDataMin(void)
{
  return this->m_RidgeOp->GetDataMin(); 
}

/**
 * Set Data Max value */
template<class TInputImage>
void
TubeExtractor2D<TInputImage>
::SetDataMax(double dataMax)
{
  this->m_RidgeOp->SetDataMin(dataMax);
  this->m_RadiusOp->SetDataMin(dataMax); 
}


/**
 * Get Data Max value */
template<class TInputImage>
double
TubeExtractor2D<TInputImage>
::GetDataMax(void)
{
  return this->m_RidgeOp->GetDataMax(); 
}

/**
 * Set Mode Retina */
template<class TInputImage>
void
TubeExtractor2D<TInputImage>
::SetModeRetina(bool modeRetina)
{
  this->m_RidgeOp->SetModeRetina(modeRetina);
  this->m_RadiusOp->SetModeRetina(modeRetina);
}


/**
 * Set Radius*/
template<class TInputImage>
void
TubeExtractor2D<TInputImage>
::SetRadius(double radius)
{
  this->m_RidgeOp->SetScale(radius);
  this->m_RadiusOp->SetRadius0(radius);
}

/**
 * Get Radius*/
template<class TInputImage>
double
TubeExtractor2D<TInputImage>
::GetRadius(void)
{
  return this->m_RidgeOp->GetScale();
}

/**
 * Extract the ridge */
template<class TInputImage>
void
TubeExtractor2D<TInputImage>
::ExtractRidge(bool extractRidge)
{
  this->m_RidgeOp->SetExtractRidge(extractRidge);
  this->m_RadiusOp->SetExtractRidge(extractRidge);
}

/**
 * Get Extract ridge */
template<class TInputImage>
bool
TubeExtractor2D<TInputImage>
::ExtractRidge(void)
{
  return this->m_RidgeOp->GetExtractRidge();
}

/**
 * Extract a valley */
template<class TInputImage>
void
TubeExtractor2D<TInputImage>
::ExtractValley(bool extractValley)
{
  this->m_RidgeOp->SetExtractValley(extractValley);
  this->m_RadiusOp->SetExtractRidge(!extractValley);
}

/**
 * Get Extract valley */
template<class TInputImage>
bool
TubeExtractor2D<TInputImage>
::ExtractValley(void)
{
  return this->m_RidgeOp->GetExtractValley();
}

/**
 * Get the ridge extractor */
template<class TInputImage>
typename RidgeExtractor2D<TInputImage>::Pointer
TubeExtractor2D<TInputImage>
::GetRidgeOp(void)
{
  return &  this->m_RidgeOp;
}
  
/**
 * Get the radius extractor */
template<class TInputImage>
typename RadiusExtractor2D<TInputImage>::Pointer
TubeExtractor2D<TInputImage>
::GetRadiusOp(void)
{
  return & this->m_RadiusOp;
}
  
/**
 * Extract the tube given the position of the first point
 * and the tube ID */
template<class TInputImage> 
bool 
TubeExtractor2D<TInputImage>
::ExtractTube(float x, float y, unsigned int tubeID)
{
 
  (*m_X)[0] = x;
  (*m_X)[1] = y;
//  this->m_RidgeOp->Print(std::cout);
  m_Tube = this->m_RidgeOp->Extract(m_X, tubeID);
//  m_Tube->Print(std::cout);
  if(!m_Tube)
    {
    if(m_Debug) std::cout << "this->m_RidgeOp->Extract() fails !" << std::endl;
    return false;
    }
    
   this->m_RidgeOp->SmoothTubeX(m_Tube, 15);
  //m_Tube->SetColor(m_Color);
   
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
  
  
  double tR = this->m_RadiusOp->GetRadius0();

  if(this->m_RidgeOp->GetAutoScale())
  {
    this->m_RadiusOp->SetRadius0(1.5*this->m_RidgeOp->GetAutoScaleUsed());
  }
  if(!this->m_RadiusOp->CalcRadii(m_Tube)) 
  {
    this->m_RadiusOp->SetRadius0(tR);
    return false;
  }
  this->m_RadiusOp->SetRadius0(tR);
      
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
TubeExtractor2D<TInputImage>
::DeleteTube(TubeType* tube)
{
  return this->m_RidgeOp->DeleteTube(tube);
}

/** 
 * Get the last tube extracted */
template<class TInputImage> 
typename TubeExtractor2D<TInputImage>::TubePointer 
TubeExtractor2D<TInputImage>
::GetLastTube(void)
{
  return m_Tube;
}
  
/**
 * Set the tube color */
template<class TInputImage> 
void
TubeExtractor2D<TInputImage>
::SetColor(float* color)
{
  m_Color = color;
}


/**
 * Get the last position */
template<class TInputImage> 
typename TubeExtractor2D<TInputImage>::VectorType* 
TubeExtractor2D<TInputImage>
::GetLastPosition(void)
{
  return m_X;
}

/**
 * Set the idle call back */
template<class TInputImage>
void
TubeExtractor2D<TInputImage>
::IdleCallBack(bool (*idleCallBack)())
{
  m_IdleCallBack = idleCallBack;
}

/**
 * Set the status callback  */
template<class TInputImage>
void
TubeExtractor2D<TInputImage>
::StatusCallBack(void (*statusCallBack)(char *, char *, int))
{
  m_StatusCallBack = statusCallBack;
  this->m_RidgeOp->StatusCallBack(statusCallBack);
  this->m_RadiusOp->SetStatusCallBack(statusCallBack);
}

/**
 * Set the status callback  */
template<class TInputImage>
void
TubeExtractor2D<TInputImage>
::NewTubeCallBack(void (*newTubeCallBack)(Tube* ))
{
  m_NewTubeCallBack = newTubeCallBack;
}

/**
 * Abort the process  */
template<class TInputImage>
void
TubeExtractor2D<TInputImage>
::AbortProcess(bool (*abortProcess)())
{
  m_AbortProcess = abortProcess;
}


}; // end namespace itk

#endif
