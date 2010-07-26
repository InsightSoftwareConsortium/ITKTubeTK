/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkRidgeExtractor.txx, v $
  Language:  C++
  Date:      $Date: 2005/10/06 16:20:06 $
  Version:   $Revision: 1.20 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRidgeExtractor_txx
#define __itkRidgeExtractor_txx

#include "itkRidgeExtractor.h"
#include "itkMatrixMath.h"
#include <itkImageRegionIterator.h>
#include <itkMinimumMaximumImageFilter.h>

namespace itk
{
 
template <class TInputImage>
class RidgeExtractorSplineValue : public UserFunc<vnl_vector<int> *, double>
{
public:

  RidgeExtractorSplineValue(RidgeExtractor<TInputImage> * newRidgeExtractor)
    {
    m_Ridge = newRidgeExtractor;
    };

  double value(vnl_vector<int> * x)
    {
    typename TInputImage::IndexType index;
    index[0] = (*x)[0];
    index[1] = (*x)[1];
    index[2] = (*x)[2];
    return m_Ridge->Intensity(index);
    };

protected:

  RidgeExtractor<TInputImage> * m_Ridge;

};


/**
 * Constructor */
template<class TInputImage>
RidgeExtractor<TInputImage>
::RidgeExtractor()
{
  m_Debug = false;
  m_Verbose = true;
   
  m_DataOp = Blur3DImageFunction<ImageType>::New();
  m_DataOp->SetScale(3); // 1.5
  m_DataOp->SetExtent(3.1); // 3

  m_StepX = 0.2;

  m_DynamicScale = false;
  m_DynamicScaleUsed = 3;
  m_RadiusExtractor = NULL;

  m_ThreshT = 0.75;
  m_ThreshX = 1.0;
  m_ThreshEV = -0.002;
  m_ThreshP2Q2 = 0.001;        // RW 0.005 // TMI 0.005 // near zero = harder
  m_ThreshP2Q2Start = 0.01;    // RW 0.01  //
  m_ThreshEVRatio = 0.2;      // RW 0.18  // TMI 0.2 // near 1 = harder
  m_ThreshEVRatioStart = 0.1;
  m_RecoveryMax = 4;

  m_SplineValueFunc = new RidgeExtractorSplineValue<TInputImage>(this);
  m_DataSpline = new SplineND(ImageDimension, m_SplineValueFunc, &m_DataSpline1D, &m_DataSplineOpt);
  m_DataSplineOpt.searchForMin(false);
  m_DataSplineOpt.tolerance(0.001);        // TMI 0.001
  m_DataSplineOpt.xStep(0.15);            // TMI 0.15
   
  m_IdleCallBack = NULL;
  m_StatusCallBack = NULL;
   
  m_Tube = NULL; 
  m_TubePointList.clear();
}

/**
 * Destructor */
template<class TInputImage>
RidgeExtractor<TInputImage>
::~RidgeExtractor()
{
  if(m_SplineValueFunc != NULL)

    {
    delete m_SplineValueFunc;
    } 
  m_SplineValueFunc = NULL;
   
  if(m_DataSpline != NULL)
    {    
    delete m_DataSpline;
    }
  m_DataSpline = NULL;

}

/**
 * Set the input image */
template<class TInputImage>
void
RidgeExtractor<TInputImage>
::SetInputImage(ImagePointer inputImage )
{
  if(m_Debug)
    {
    std::cout << std::endl << "Ridge::SetInputImage" << std::endl;
    }

  m_Image = inputImage;

  if(m_Image) 
    {
    m_DataOp->SetInputImage(inputImage);
    
    typedef MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
    typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
    minMaxFilter->SetInput(m_Image);
    minMaxFilter->Update();
    m_DataMin = minMaxFilter->GetMinimum();
    m_DataMax = minMaxFilter->GetMaximum();
    m_DataRange = m_DataMax-m_DataMin;

    if(m_Debug)
      {
      std::cout << "  Minimum = " << m_DataMin << std::endl;
      std::cout << "  Maximum = " << m_DataMax << std::endl;
      std::cout << "  Data Range = " << m_DataRange << std::endl;
      }

    unsigned int size[ImageDimension];
    for(unsigned int i=0; i<ImageDimension; i++)
      {
      size[i]= m_Image->GetLargestPossibleRegion().GetSize()[i];
      m_ExtractBoundMax[i] = size[i]-1;
      m_ExtractBoundMin[i] = (int)0.0;
      }

    vnl_vector<int> v(3);
    v[0] = m_ExtractBoundMin[0];
    v[1] = m_ExtractBoundMin[1];
    v[2] = m_ExtractBoundMin[2];
    m_DataSpline->xMin(v);
    v[0] = m_ExtractBoundMax[0];
    v[1] = m_ExtractBoundMax[1];
    v[2] = m_ExtractBoundMax[2];
    m_DataSpline->xMax(v);

    /** Allocate the mask image */
    m_DataMask = MaskType::New();

    typename ImageType::RegionType region;
    m_DataMask->SetRegions( m_Image->GetLargestPossibleRegion() );
    m_DataMask->CopyInformation( m_Image );
    m_DataMask->Allocate();
    m_DataMask->FillBuffer(0);

    } // end Image == NULL
}

/**
 * Set Data Min value */
template<class TInputImage>
void
RidgeExtractor<TInputImage>
::SetDataMin(double dataMin)
{
  m_DataMin = dataMin;
  m_DataRange = m_DataMax-m_DataMin;
}

/**
 * Set Data Min value */
template<class TInputImage>
void
RidgeExtractor<TInputImage>
::SetDataMax(double dataMax)
{
  m_DataMax = dataMax;
  m_DataRange = m_DataMax-m_DataMin;
}

/**
 * Set the scale */
template<class TInputImage>
void
RidgeExtractor<TInputImage>
::SetScale(double scale)
{
  if(m_Verbose)
    {
    std::cout << "Ridge::SetScale = " << scale << std::endl;
    }
  m_DataSpline->newData(true);
  m_DataOp->SetScale(scale);
}

/**
 * Get the scale */
template<class TInputImage>
double
RidgeExtractor<TInputImage>
::GetScale(void)
{
  return m_DataOp->GetScale();
}

/**
 * Set the extent */
template<class TInputImage>
void
RidgeExtractor<TInputImage>
::SetExtent(double extent)
{
  m_DataSpline->newData(true);
  m_DataOp->SetExtent(extent);
}


/**
 * Get the extent */
template<class TInputImage>
double
RidgeExtractor<TInputImage>
::GetExtent(void)
{
  return m_DataOp->GetExtent();
}

/**
 * Get the data spline */
template<class TInputImage>
SplineND* 
RidgeExtractor<TInputImage>   
::GetDataSpline(void)
{
  return m_DataSpline;
}
  
/**
 * Get the data spline 1D*/
template<class TInputImage>
Spline1D* 
RidgeExtractor<TInputImage>
::GetDataSpline1D(void)
{
  return & m_DataSpline1D;
}
  
/**
 * Get the data spline optimizer*/
template<class TInputImage>
Optimizer1D* 
RidgeExtractor<TInputImage>
::GetDataSplineOptimizer(void)
{
  return & m_DataSplineOpt;
}


/**
 * Set the dynamic scale */
template<class TInputImage>
void
RidgeExtractor<TInputImage>
::SetDynamicScale(bool dynamicScale)
{
  if(m_RadiusExtractor)
    {
    this->m_DynamicScale = dynamicScale;
    }
}


/**
 * Set the radius extractor */
template<class TInputImage>
void 
RidgeExtractor<TInputImage>
::SetRadiusExtractor(RadiusExtractor<TInputImage> * radiusExtractor)
{
  m_RadiusExtractor = radiusExtractor;
}

/**
 * Return the intensity */
template<class TInputImage>
double 
RidgeExtractor<TInputImage>
::Intensity(IndexType & x)
{
  double tf = (m_DataOp->EvaluateAtIndex(x)-m_DataMin)/m_DataRange;
  
  if(tf<0)
    {
    tf = 0;
    }
  else if(tf>1)
    {
    tf = 1;
    }
  
  return tf;
}

/**
 * Traverse one way
 */
template<class TInputImage>
double 
RidgeExtractor<TInputImage>
::Ridgeness(ContinuousIndexType & x)
{    
  if(m_Debug)
    {
    std::cout << "Ridge::Ridgeness" << std::endl;
    }

  m_X = x;
  m_XVal = m_DataSpline->valueJet(m_X.GetVnlVector(),
                                  m_XD.GetVnlVector(),
                                  m_XH.GetVnlMatrix().as_ref());
   
  // HACK
  unsigned int flat = 0;
  for(unsigned int i=0;i<ImageDimension;i++)
    {
    if(fabs(m_XD[i]) < 1e-10) 
      {
      flat++;
      }
    }

  if(flat == ImageDimension)
    {
    if(m_Verbose)
      {
      std::cout << "Ridgeness: Warning: Encountered flat region" << std::endl;
      }
    m_XP = m_XQ = m_XR = 1;
    // HACK
    return 2;
    }

  Eigen(m_XH.GetVnlMatrix().as_ref(),
        m_XHEVect.GetVnlMatrix().as_ref(), m_XHEVal.GetVnlVector(), false);

  m_XP = m_XHEVect(0, 0)*m_XD[0] 
         + m_XHEVect(1, 0)*m_XD[1] 
         + m_XHEVect(2, 0)*m_XD[2];
  m_XQ = m_XHEVect(0, 1)*m_XD[0] 
         + m_XHEVect(1, 1)*m_XD[1] 
         + m_XHEVect(2, 1)*m_XD[2];
  m_XR = m_XHEVect(0, 2)*m_XD[0] 
         + m_XHEVect(1, 2)*m_XD[1] 
         + m_XHEVect(2, 2)*m_XD[2];

  m_XD.Normalize();

  if(m_Debug) 
    {
    std::cout << "  m_X = " << m_X[0] << ", " << m_X[1] << ", " << m_X[2] << std::endl;
    std::cout << "  Scale = " << m_DataOp->GetScale() << std::endl;
    std::cout << "  m_XVal = " << m_XVal << std::endl;
    std::cout << "  m_XP = " << m_XP << std::endl;
    std::cout << "  m_XQ = " << m_XQ << std::endl;
    std::cout << "  m_XR = " << m_XR << std::endl;
    std::cout << "  m_XD = " 
              << m_XD[0] << ", " << m_XD[1] << ", " << m_XD[2] << std::endl;
    std::cout << "  m_XHEVal = " 
              << m_XHEVal[0] << ", " << m_XHEVal[1] << ", " << m_XHEVal[2] 
              << std::endl;
    std::cout << "  m_XHEVect(0) = " 
              << m_XHEVect(0, 0) << ", " << m_XHEVect(1, 0) << ", " 
              << m_XHEVect(2, 0) << std::endl;
    std::cout << "  m_XHEVect(1) = " << m_XHEVect(0, 1) << ", " 
              << m_XHEVect(1, 1) << ", " << m_XHEVect(2, 1) << std::endl;
    std::cout << "  m_XHEVect(2) = " << m_XHEVect(0, 2) << ", " 
              << m_XHEVect(1, 2) << ", " << m_XHEVect(2, 2) << std::endl;
    std::cout << "  m_XHEVect(0)*m_XHEVect(1) = " 
              << m_XHEVect(0, 0)*m_XHEVect(0, 1) +
                 m_XHEVect(1, 0)*m_XHEVect(1, 1) +
                 m_XHEVect(2, 0)*m_XHEVect(2, 1) 
             << std::endl;
    std::cout << "  Result = p^2 + q^2 = " << m_XP*m_XP+m_XQ*m_XQ << std::endl;
    }   
  

  return m_XP*m_XP+m_XQ*m_XQ;
}

/**
 * Traverse one way
 * Need to be implemented in 2D */
template<class TInputImage>
typename RidgeExtractor<TInputImage>::TubeType * 
RidgeExtractor<TInputImage>
::TraverseOneWay(ContinuousIndexType & newX, VectorType & newT,
                 NormalPlaneMatrixType & newN, int dir)
{
  if(m_Debug)
    {
    std::cout << "Ridge::TraverseOneWay" << std::endl;
    }

  double                 lVal;
  ContinuousIndexType    lX;
  VectorType             lT;
  NormalPlaneMatrixType  lN;
  VectorType             lNTEVal;
  VectorType             lStepDir;
  NormalPlaneMatrixType  lSearchDir;
  
  ContinuousIndexType    pX;
  VectorType             pT;
  NormalPlaneMatrixType  pN;
  VectorType             pStepDir;
  NormalPlaneMatrixType  pSearchDir;
  
  lX = newX;
  lT = newT;
  lN = newN;
  lStepDir = newT;
  lSearchDir = newN;
  
  pX = lX;
  pT = lT;
  pN = lN;
  pStepDir = lStepDir;
  pSearchDir = lSearchDir;
  unsigned int i, j, k, l;
  
  VectorType prod;
  VectorType tV;
  
  typename ImageType::IndexType index;
  for(i=0; i<ImageDimension; i++)
    {
    index[i] = (int)(lX[i]+0.5);
    if(index[i]<(m_ExtractBoundMin)[i] ||
       index[i]>(m_ExtractBoundMax)[i])
      {
      if(m_Verbose)
        {
        std::cout << "Ridge: TraverseOneWay: Exited boundary" << std::endl;
        }
      return m_Tube;
      }
    }

  double p2q2;
  
  typename MaskType::PixelType value = m_DataMask->GetPixel(index);
  if(value != 0 && (int)value != m_TubeID)
    {
    if(m_Verbose)
      {
      std::cout << "Ridge: TraverseOneWay: Encountered another tube" << std::endl;
      }
    return m_Tube;
    }
  else
    {
    m_DataMask->SetPixel(index, (PixelType)(m_TubeID+(m_TubePointCount/10000.0)));
    if(dir == 1)
      {
      if(m_Debug)
        {
        std::cout << "Ridge: dir = 1" << std::endl;
        }
      // p2q2 = Ridgeness(&(lX));
      TubePointType pnt;
      pnt.SetID(m_TubePointCount);
      pnt.SetPosition(lX[0], lX[1], lX[2]);
      pnt.SetNormal1(m_XHEVect(0,0), m_XHEVect(1,0), m_XHEVect(2,0));
      pnt.SetNormal2(m_XHEVect(0,1), m_XHEVect(1,1), m_XHEVect(2,1));
      pnt.SetAlpha1(m_XHEVal[0]);
      pnt.SetAlpha2(m_XHEVal[1]);
      if(m_XHEVal[0]!=0)
        {
        pnt.SetRidgeness(m_XHEVal[1]/m_XHEVal[0]);
        }
      else
        {
        pnt.SetRidgeness(0.0);
        }
      //m_Tube->GetPoints().push_back(pnt);
      m_TubePointList.push_front(pnt);
      m_TubePointCount++;
      }
    }
     
  double iScale0 = GetScale();
  
  int pSize = m_Tube->GetPoints().size();
  
  double stepFactor0 = 1;
  double stepFactor = stepFactor0;
  
  double tf;
  int recovery = 0;
  int prevRecoveryPoint = -(1.5/m_StepX);
  while(recovery<m_RecoveryMax 
        && prevRecoveryPoint+(2/m_StepX)>m_TubePointCount) 
    {
    if(recovery > 0) 
      {
      if(m_Verbose)
        {
        std::cout << "Attempting recovery : " << recovery;
        std::cout << " : Scale = " << GetScale() << std::endl;
        std::cout << "   x=" << pX[0] << ", " << pX[1] << ", " << pX[2] << std::endl;
        }
      switch(recovery)
        {
        default:
        case 1:
          prevRecoveryPoint = m_TubePointCount;
          stepFactor = 1.5 * stepFactor0;
          break;
        case 2:
          stepFactor = 0.5 * stepFactor0;
          SetScale(GetScale() * 1.25);
          break;
        case 3:
          stepFactor = 2.0 * stepFactor0;
          SetScale(GetScale() * 1.25);
          break;
        }
      if(m_Verbose)
        {
        std::cout << "   Point = " << m_TubePointCount 
                  << ": Recovery: new scale = " << GetScale() << std::endl;
        }    
      lX = pX;
      lT = pT;
      lN = pN;
      lStepDir = pStepDir;
      lSearchDir = pSearchDir;
      }
    else 
      {
      if(m_Debug)
        {
        std::cout << "Ridge: No recovery needed" << std::endl;
        }
      if(!m_DynamicScale && GetScale() != iScale0)
        {
        SetScale(iScale0+0.5*(iScale0-GetScale()));
        }

      if(fabs(dot_product(lStepDir.GetVnlVector(), pStepDir.GetVnlVector()))<1-0.5*(1-m_ThreshT)) 
        {
        if(m_Verbose)
          {
          std::cout << "Curving" << std::endl;
          }  
        stepFactor *= (double)0.75;
        }
      else
        {   
        if(stepFactor<stepFactor0)
          {
          stepFactor *= 1.25;
          }
        } 
      if(stepFactor>stepFactor0)
        {   
        stepFactor = stepFactor0;
        }
      pX = lX;
      pT = lT;
      pN = lN;
      pStepDir = lStepDir;
      pSearchDir = lSearchDir;
      }
    
    vnl_vector<double> v = ComputeLineStep(lX.GetVnlVector(), m_StepX*stepFactor, lStepDir.GetVnlVector());
    lX[0] = v[0];
    lX[1] = v[1];
    lX[2] = v[2];
    if(m_Debug)
      {
      std::cout << "Ridge: Computed line step = " << v << std::endl;
      std::cout << "Ridge: TraverseOW: lStepDir = " 
                << lStepDir[0] << ", " << lStepDir[1] << ", " 
                << lStepDir[2] << std::endl;
      std::cout << "Ridge: TraverseOW: lX0 = " 
                << lX[0] << ", " << lX[1] << ", " << lX[2] << std::endl;
      std::cout << "Ridge: TraverseOW: lSearchDir1 = " 
                << lSearchDir[0][0] << ", " << lSearchDir[1][0] << ", " 
                << lSearchDir[2][0] << std::endl;
      std::cout << "Ridge: TraverseOW: lSearchDir2 = " 
                << lSearchDir[0][1] << ", " << lSearchDir[1][1] << ", " 
                << lSearchDir[2][1] << std::endl;
      }
    
    if(!m_DataSpline->extreme(lX.GetVnlVector(), &lVal, 2, lSearchDir.GetVnlMatrix().as_ref()))
      {
      if(m_Verbose)
        {
        std::cout << "*** Ridge terminated: Local max not found" << std::endl;
        }
      recovery++;
      continue;
      }

    if(m_Debug)
      {
      std::cout << "Ridge: TraverseOW: lXExtreme = " 
                << lX[0] << ", " << lX[1] << ", " << lX[2] << std::endl;
      }

    if((lX[0])<(m_ExtractBoundMin)[0] ||
       (lX[0]+0.5)>(m_ExtractBoundMax)[0] ||
       (lX[1])<(m_ExtractBoundMin)[1] ||
       (lX[1]+0.5)>(m_ExtractBoundMax)[1] ||
       (lX[2])<(m_ExtractBoundMin)[2] ||
       (lX[2]+0.5)>(m_ExtractBoundMax)[2]) 
      {
      if(m_Verbose)
        {
        std::cout << "*** Ridge term: Exited extraction bounds" << std::endl;
        }
      break;
      }
    
    p2q2 = Ridgeness(lX);  
    
    for(i=0; i<3; i++) 
      {
      for(j=0; j<3; j++)
        {
        tV[j] = m_XHEVect(j, i);
        }
      prod[i] = fabs(dot_product(lStepDir.GetVnlVector(), tV.GetVnlVector()));
      }

    k = 0;
    tf = prod[0];
    for(i=1; i<3; i++) 
      {
      if(prod[i]>tf) 
        {
        k = i;
        tf = prod[i];
        }
      }
    if(k != 2)
      {
      if(m_Verbose)
        {
        std::cout << "Mixing things up: Chosen t=evect#" 
                  << k << " tf=" << tf << std::endl;
        }
      }
    j = (k+1)%3;
    l = (k+2)%3;
    if(m_XHEVal[j]>m_XHEVal[l]) 
      {
      i = j;
      j = l;
      l = i;
      }

    for(i=0; i<3; i++)
      {
      lN[i][0] = m_XHEVect(i, j);
      }   
    lNTEVal[0] = m_XHEVal[j];
         
    for(i=0; i<3; i++)
      {
      lN[i][1] = m_XHEVect(i, l);
      }  
    lNTEVal[1] = m_XHEVal[l];
    lSearchDir = lN;
         
    for(i=0; i<3; i++)
      {
      lT[i] = m_XHEVect(i, k);
      } 
    lNTEVal[2] = m_XHEVal[k];
         
    tf = dot_product(pT.GetVnlVector(), lT.GetVnlVector());
    if(tf<0) 
      {
      lT[0] = lT[0] * -1;
      lT[1] = lT[1] * -1;
      lT[2] = lT[2] * -1;
      }
         
    lStepDir[0] = (pStepDir[0] * 0.25) + (lT[0] * 0.75);
    lStepDir[1] = (pStepDir[1] * 0.25) + (lT[1] * 0.75);
    lStepDir[2] = (pStepDir[2] * 0.25) + (lT[2] * 0.75);

    tf = dot_product(lStepDir.GetVnlVector(), pStepDir.GetVnlVector());
    if(m_ThreshT>0 && fabs(tf)<m_ThreshT) 
      {
      if(m_Verbose)
        {
        std::cout << "*** Ridge term: Rapid change in step direction "
                  << "(" << fabs(tf) << ")" << std::endl;
        }
      recovery++;
      continue;
      }
     
    tf = sqrt(ComputeEuclideanDistanceVector(lX.GetVnlVector(), pX.GetVnlVector()));
    if(m_ThreshX>0 && tf>m_ThreshX + 0.1*recovery) 
      {
      if(m_Verbose)
        {
        std::cout << "*** Ridge term: Rapid change in spatial location "
                  << "(" << tf << ")" << std::endl;
        } 
      recovery++;
      continue;
      }
       
    if(p2q2>m_ThreshP2Q2) 
      {
      if(m_Verbose)
        {
        std::cout << "*** Ridge terminated: Local max not a ridge point " 
                  << "(p2q2 = " << p2q2 << ")" << std::endl;
        }
      tf = lNTEVal[1]/lNTEVal[0];
         
      if(m_Debug)
        {
        std::cout << "       EV Ratio = " << tf << std::endl;
        }
      recovery++;
      continue;
      }

    if(lNTEVal[1]>m_ThreshEV)
      {
      if(m_Verbose)
        {
        std::cout << "*** Ridge terminated: EV thresh: low curvature "
                  << "(" << lNTEVal[0] << ", " << lNTEVal[1] << ")" << std::endl;
        std::cout << "       P2Q2 = " << p2q2 << std::endl;
        }
      if(fabs(lNTEVal[0])!=0 && p2q2 != 0)
        {
        recovery++;
        }
      else
        {       
        recovery = m_RecoveryMax;
        }   
      continue;
      }

    tf = lNTEVal[1]/lNTEVal[0];
    if(tf<m_ThreshEVRatio) 
      {
      if(m_Verbose)
        {
        std::cout << "*** Ridge terminated: EV Ratio : Planar point " 
                  << "(" << tf << ")" << std::endl;
        std::cout << "       P2Q2 = " << p2q2 << std::endl;
        }
      if(fabs(lNTEVal[0])!=0 && p2q2 != 0)
        {
        recovery++;
        }   
      else
        {
        recovery = m_RecoveryMax;
        }   
      continue;
      }
 
    for(i=0; i<ImageDimension; i++)
      {
      index[i] = (int)(lX[i]+0.5);
      }
    tf = m_DataMask->GetPixel(index);
   
    if(tf != 0) 
      {
      if((int)tf != m_TubeID ||
         m_TubePointCount-((tf-(int)tf)*10000)>20/m_StepX) 
        {
        if(m_Verbose)
          {
          std::cout << "*** Ridge terminated: Revisited self" << std::endl;
          }
        break;
        }
      }
    else
      {
      m_DataMask->SetPixel(index, (PixelType)(m_TubeID+(m_TubePointCount/10000.0)));
      }   

    /** Show the satus every 50 points */
    if(m_TubePointCount%50==0)
      {
      if(m_StatusCallBack) 
        {
        char st[80];
        sprintf(st, "Point #%d", m_TubePointCount);
        m_StatusCallBack(NULL, st, 0);
        }
      }      
    if(m_Debug)
      {
      std::cout << "Ridge: TraverseOW: Adding point " << m_TubePointCount << " = " 
                << lX[0] << ", " << lX[1] << ", " << lX[2] 
                << std::endl;
      }

    TubePointType pnt;
    pnt.SetID(m_TubePointCount);
    pnt.SetPosition(lX[0], lX[1], lX[2]);
    pnt.SetNormal1(lN(0,0), lN(1,0), lN(2,0));
    pnt.SetAlpha1(lNTEVal[0]);
    pnt.SetNormal2(lN(0,1), lN(1,1), lN(2,1));
    pnt.SetAlpha2(lNTEVal[1]);

    if(lNTEVal[0] !=0)
      {
      pnt.SetRidgeness(lNTEVal[1]/lNTEVal[0]);
      }
    else
      {
      pnt.SetRidgeness(0.0);
      }

    if(dir == 1)
      {   
      m_TubePointList.push_front(pnt);
      }   
    else
      {   
      m_TubePointList.push_back(pnt);
      } 

    m_TubePointCount++;

    recovery = 0;
    prevRecoveryPoint = m_TubePointCount;

    if(m_TubePointCount/25.0 == m_TubePointCount/25) 
      {
      if(m_IdleCallBack)
        {
        m_IdleCallBack();
        }    
      if(m_DynamicScale && m_TubePointCount/50.0 == m_TubePointCount/50)
        {
        if(m_Debug)
          {
          std::cout << "Ridge: TraverseOW: DynamicScale" << std::endl;
          }
        if(m_RadiusExtractor) 
          {
          TubePointType tmpPoint;
          pX[0] = (pX[0]+lX[0])/2;
          pX[1] = (pX[1]+lX[1])/2;
          pX[2] = (pX[2]+lX[2])/2;
          tmpPoint.SetPosition(pX[0], pX[1], pX[2]);
          tmpPoint.SetTangent(lStepDir[0], lStepDir[1], lStepDir[2]);
          tmpPoint.SetRadius(m_DynamicScaleUsed);
          m_RadiusExtractor->SetRadius0(m_DynamicScaleUsed);
          if(m_RadiusExtractor->CalcOptimalScale(tmpPoint))
            {
            m_DynamicScaleUsed = (2*tmpPoint.GetRadius()+m_DynamicScaleUsed)/3;// use to avg with cDy)/2.0;
            }
          if(m_DynamicScaleUsed<0.5)
            {
            m_DynamicScaleUsed = 0.5;
            }
          if(m_StatusCallBack) 
            {
            char s[80];
            sprintf(s, "Extract: Ridge: DS = %1.1f", 
              m_DynamicScaleUsed);
            m_StatusCallBack(s, NULL, 0);
            }
          else if(m_Verbose)
            {   
            std::cout << "Dynamic Scale = " << m_DynamicScaleUsed << std::endl;
            }
          SetScale(m_DynamicScaleUsed);
          m_RadiusExtractor->SetRadius0(m_DynamicScaleUsed);
          }
        }
      } 
    } 
  
  if(m_Verbose)
    {
    std::cout << "*** Ridge terminated: Cannot recover" << std::endl;
    std::cout << "    Length = " << m_Tube->GetPoints().size()-pSize << std::endl;
    }
  SetScale(iScale0);

  return m_Tube;
}
  
/**
 * Compute the local ridge
 */
template<class TInputImage>
bool
RidgeExtractor<TInputImage>
::LocalRidge(ContinuousIndexType & newX)
{  
  if(m_Debug)
    {
    std::cout << "Ridge::LocalRidge" << std::endl;
    }

  typename ImageType::IndexType index;
  for(unsigned int i=0;i<ImageDimension;i++)
    {
    index[i] = (int)(newX[i]);
    if(index[i]<0 
       || index[i]>=m_Image->GetLargestPossibleRegion().GetSize()[i])
      {
      if(m_StatusCallBack)
        { 
        m_StatusCallBack(NULL, "Exited Image", 0);
        }
      if(m_Verbose)
        {
        std::cout << "RidgeExtractor::LocalRidge() : Exited Image 2" << std::endl;
        }
      return false;
      }
    }

  if(Ridgeness(newX)>=2)
    { 
    if(m_StatusCallBack) 
      {
      m_StatusCallBack(NULL, "RidgeExtractor::LocalRidge(): Ridgeness>=2", 0);
      }
    if(m_Verbose)
      {
      std::cout << "RidgeExtractor::LocalRidge() : Ridgeness>=2" << std::endl;
      }
    return false;
    }

  double                tf;  
  double                val;
  NormalPlaneMatrixType lN;
  ContinuousIndexType   pX;
   
  int loop;
  double p2q2 = 0;
  // Gradient Ascent
  for(loop=0; loop<ImageDimension-1; loop++)
    {
    pX = newX;
    for(int i=0; i<ImageDimension; i++)
      {  
      lN(i, 0) = m_XD[i];
      }
  
    m_DataSpline->extreme(newX.GetVnlVector(), &val, 1, lN.GetVnlMatrix().as_ref());
    
    for(int i=0; i<ImageDimension; i++)
      {
      index[i] = (int)(newX[i]);
      if(index[i]<0 
         || index[i]>=m_Image->GetLargestPossibleRegion().GetSize()[i])
        {
        if(m_StatusCallBack)
          {
          m_StatusCallBack(NULL, "Exited Image", 0);
          }
        if(m_Verbose)
          {
          std::cout << "RidgeExtractor::LocalRidge() : Exited Image 3" << std::endl;
          }
        return false;
        }
      }

    if(Ridgeness(newX)>=2)
      { 
      if(m_StatusCallBack) 
        {
        m_StatusCallBack(NULL, "RidgeExtractor::LocalRidge(): Ridgeness>=2", 0);
        }
      if(m_Verbose)
        {
        std::cout << "RidgeExtractor::LocalRidge() : Ridgeness>=2" << std::endl;
        }
      return false;
      }

    if(m_DataMask->GetPixel(index) != 0)
      {
      if(m_StatusCallBack) 
        {
        m_StatusCallBack(NULL, "Revisited voxel", 0);
        }
      if(m_Verbose)
        {
        std::cout << "RidgeExtractor::LocalRidge() : Revisited voxel 1" << std::endl;
        }
      return false;
      }
    }

  for(loop=0; loop<ImageDimension; loop++)
    {
    for(int i=0; i<ImageDimension; i++) 
      {
      lN(i, 0) = m_XHEVect(i, 0);
      lN(i, 1) = m_XHEVect(i, 1);
      }

    // Local 1D Ridge
    m_DataSpline->extreme(newX.GetVnlVector(), &val, 2, lN.GetVnlMatrix().as_ref());

    for(unsigned int i=0;i<ImageDimension;i++)
      {
      index[i]=(int)(newX[i]);
      if(index[i]<0 || index[i]>=m_Image->GetLargestPossibleRegion().GetSize()[i])
        {
        if(m_StatusCallBack)
          {
          m_StatusCallBack(NULL, "Exited Image", 0);
          }
        if(m_Verbose) 
          {
          std::cout << "RidgeExtractor::LocalRidge() : Exited Image 5" << std::endl;
          }
        return false;
        }
      }

    if(m_DataMask->GetPixel(index) != 0)
      {
      if(m_StatusCallBack)
        {
        m_StatusCallBack(NULL, "Revisited voxel", 0);
        }
      if(m_Verbose)
        {
        std::cout << "RidgeExtractor::LocalRidge() : Revisited voxel 3" << std::endl;
        }
      return false;
      } 

    p2q2 = Ridgeness(newX);
    if(p2q2 > m_ThreshP2Q2Start) 
      {   
      continue;
      }

    tf = m_XHEVal[1]/m_XHEVal[0];
    if(tf < m_ThreshEVRatioStart)
      {
      continue;
      }

    return true;
    }


  if(p2q2 > m_ThreshP2Q2Start) 
    {
    if(m_StatusCallBack) 
      {
      m_StatusCallBack(NULL, "p2q2 failure", 0); 
      }
    if(m_Verbose) 
      {
      std::cout << "LocalRidge : p2q2 failure" << std::endl; 
      }
    }

  if(tf < m_ThreshEVRatioStart)
    {   
    if(m_StatusCallBack) 
      {
      m_StatusCallBack(NULL, "EVRatio failure", 0); 
      }
    if(m_Verbose) 
      {
      std::cout << "LocalRidge : EVRatio failure" << std::endl; 
      }
    }

  return false;
}

/**
 * Extract a tube 
 */
template<class TInputImage>
typename RidgeExtractor<TInputImage>::TubePointer
RidgeExtractor<TInputImage>
::Extract(ContinuousIndexType & newX, int tubeID)
{
  ContinuousIndexType lX;
  lX = newX;

  double scaleOriginal = GetScale();
  double scale0 = scaleOriginal;
  double radiusOriginal = m_RadiusExtractor->GetRadius0();

  if(!LocalRidge(lX))
    {
    if(m_Verbose) 
      {
      std::cout << "LocalRidge fails at " 
                << lX[0] << " : " << lX[1] << " : " << lX[2] << std::endl;
      }
    return NULL;
    }
  
  if(m_Debug)
    {
    std::cout << "*** Ridge found at " 
              << lX[0] << ", " << lX[1] << ", " << lX[2] << std::endl;
    }

  m_TubePointList.clear();

  NormalPlaneMatrixType lN;
  VectorType lT;
   
  if(m_DynamicScale && m_RadiusExtractor != NULL) 
    {
    for(int i=0; i<ImageDimension; i++)
      {
      lT[i] = m_XHEVect(i, ImageDimension-1);
      } 
    TubePointType tmpPoint;
    tmpPoint.SetPosition(lX[0], lX[1], lX[2]);
    tmpPoint.SetTangent(lT[0], lT[1], lT[2]);
    tmpPoint.SetRadius(scale0);
    if(!m_RadiusExtractor->CalcOptimalScale(tmpPoint, true)) 
      {
      if(m_Debug && m_StatusCallBack)
        {
        m_StatusCallBack("Extract: Ridge: AS = ?", "Error: Medial Max Not Found", 0);
        }
      m_DynamicScaleUsed = scale0;
      }
    else 
      {
      m_DynamicScaleUsed = (2*tmpPoint.GetRadius()+GetScale())/3;
      }
    if(m_DynamicScaleUsed<0.5)
      {
      m_DynamicScaleUsed = 0.5;
      }

    SetScale(m_DynamicScaleUsed);
    m_RadiusExtractor->SetRadius0(m_DynamicScaleUsed);
    if(m_Debug)
      {   
      std::cout << "DynamicScale = " << GetScale() << std::endl;
      }
    for(int i=0; i<ImageDimension; i++)
      {
      lX[i] = (lX[i] + newX[i])/2;
      }
    if(!LocalRidge(lX)) 
      {
      if(m_StatusCallBack)
        {
        m_StatusCallBack("AS Failure", NULL, 0);
        }
      if(m_Verbose)
        {
        std::cout << "RidgeExtractor:Extract(): AS Failure" << std::endl;
        }
      m_DynamicScaleUsed = scaleOriginal;
      SetScale(scaleOriginal);
      m_RadiusExtractor->SetRadius0(radiusOriginal);
      return NULL;
      }
    scale0 = m_DynamicScaleUsed;
    SetScale(scale0);
    m_RadiusExtractor->SetRadius0(scale0);
    }  

  m_Tube = TubeType::New();
  m_TubeID = tubeID;
  m_Tube->SetId(tubeID);
  m_TubePointCount = 0;
  
  for(int i=0; i<ImageDimension; i++) 
    {
    lN[i][0] = m_XHEVect(i, 0);
    lN[i][1] = m_XHEVect(i, 1);
    lT[i] = m_XHEVect(i, ImageDimension-1);
    }
  TraverseOneWay(lX, lT, lN, 1);

  if(m_DynamicScale)
    {
    SetScale(scale0);
    m_RadiusExtractor->SetRadius0(scale0);
    }
  
  for(unsigned int i=0;i<ImageDimension;i++)
    {
    lT[i] = -1 * lT[i];
    }

  TraverseOneWay(lX, lT, lN, -1);

  // Set the list of tubepoints
  std::list<TubePointType>::const_iterator ptIt = m_TubePointList.begin();
  
  while( ptIt != m_TubePointList.end())
    {
    m_Tube->GetPoints().push_back(*ptIt);
    ptIt++;
    }

  // return to user defaults
  if(m_DynamicScale)
    {
    SetScale(scaleOriginal);
    m_RadiusExtractor->SetRadius0(radiusOriginal);
    }
     
  if(m_Tube->GetPoints().size()<10) 
    {
    if(m_StatusCallBack)
      {
      m_StatusCallBack("Extract: Ridge", "Too short", 0);
      }
    DeleteTube(m_Tube);
    m_Tube = NULL;
    return m_Tube;
    }
   
  if(m_Verbose)
    {
    std::cout << "*** Extracted ridge of " << m_Tube->GetPoints().size() 
              << " points." << std::endl;
    }
  
  //
  // Calculate tangents
  //
  if(m_Tube && m_Tube->GetPoints().size() > 0)
    {
    VectorType tangent;
    tangent.Fill(0.0);
    tangent[0] = 1;
  
    std::vector<TubePointType>::iterator i, j, k;
    i = m_Tube->GetPoints().begin();
    k = m_Tube->GetPoints().end();
    k--;
    while(i != k)
      {
      j = i;
      ++j;
      tangent = (*j).GetPosition() - (*i).GetPosition();
      (*i).SetTangent(tangent);
      ++i;
      }
     
    (*k).SetTangent(tangent);
    } 

  if(m_StatusCallBack) 
    {
    char s[80];
    sprintf(s, "%d points", m_Tube->GetPoints().size());
    m_StatusCallBack("Extract: Ridge", s, 0);
    }
  return m_Tube;
}


/**
 * Smooth a tube */
template<class TInputImage>
void
RidgeExtractor<TInputImage>
::SmoothTubeX(TubeType * tube, int h)
{
  TubeType::PointType avg;

  //std::vector<TubePointType> & points = tube->GetPoints();
  //std::vector<TubePointType>::iterator pnt, pntT;
  
  TubeType::PointListType &points = tube->GetPoints();
  TubeType::PointListType::iterator pnt, pntT;
  
  //std::vector<TubePointType>::iterator begin = points.begin();
  //std::vector<TubePointType>::iterator end = points.end();

  TubeType::PointListType::iterator begin = points.begin();
  TubeType::PointListType::iterator end = points.end();

  DeleteTube(tube);

  for(pnt = begin; pnt != end; pnt++)
    {
    int cnt = 0;
    avg.Fill(0);

    if(pnt != begin)
      {
      pntT = pnt;
      for(int i=0; i<h/2 && pntT!=begin; i++, pntT--)
        {
        avg[0] += (*pntT).GetPosition()[0];
        avg[1] += (*pntT).GetPosition()[1];
        avg[2] += (*pntT).GetPosition()[2];
        cnt++;
        }
      }
    if(pnt != end)
      {
      pntT = pnt;
      pntT++;
      for(int i=0; i<h/2 && pntT!=tube->GetPoints().end(); i++, pntT++)
        {
        avg[0] += (*pntT).GetPosition()[0];
        avg[1] += (*pntT).GetPosition()[1];
        avg[2] += (*pntT).GetPosition()[2];
        cnt++;
        }
      }
    if(cnt>0)
      {
      for(unsigned int i=0;i<ImageDimension;i++)
        {
        avg[i] /= cnt;
        }
      (*pnt).SetPosition(avg);
      }
    }

  AddTube(tube);
}


/** Delete a tube 
 *  FIX ME: Need to be reimplemented in 2D */
template<class TInputImage>
bool
RidgeExtractor<TInputImage>
::DeleteTube(TubeType * tube)
{    
  int x, y, z;
  int i, j, k;
  double r;
  std::vector<TubePointType>::iterator pnt;
  for(pnt = tube->GetPoints().begin(); pnt != tube->GetPoints().end(); pnt++)
    {
    x = (int)((*pnt).GetPosition()[0]+0.5);
    y = (int)((*pnt).GetPosition()[1]+0.5);
    z = (int)((*pnt).GetPosition()[2]+0.5);
    
    if(x>=0 && x<m_DataMask->GetLargestPossibleRegion().GetSize()[0] &&
       y>=0 && y<m_DataMask->GetLargestPossibleRegion().GetSize()[1] &&
       z>=0 && z<m_DataMask->GetLargestPossibleRegion().GetSize()[2])
      {
      typename ImageType::IndexType index;
      index[0]=x;
      index[1]=y;
      index[2]=z;
      m_DataMask->SetPixel(index, 0);
      r = (*pnt).GetRadius();
      r *= 1.25;
      if(r>1 &&
         x-r>=0 && x+r<m_DataMask->GetLargestPossibleRegion().GetSize()[0] &&
         y-r>=0 && y+r<m_DataMask->GetLargestPossibleRegion().GetSize()[1] &&
         z-r>=0 && z+r<m_DataMask->GetLargestPossibleRegion().GetSize()[2])
        {
        for(i=-(int)r; i<=(int)r; i++)
          {
          for(j=-(int)r; j<=(int)r; j++)
            {
            for(k=-(int)r; k<=(int)r; k++)
              {
              index[0]=x+k;
              index[1]=y+j;
              index[2]=z+i;
              m_DataMask->SetPixel(index, 0);
              }
            }
          }
        }
      } 
    }
  return true;
}

/**
 * Add a tube */
template<class TInputImage>
bool
RidgeExtractor<TInputImage>
::AddTube(TubeType * tube)
{
  m_TubeID = tube->GetId();
  m_TubePointCount = 0;
  int x, y, z, i, j, k;
  double r;
  std::vector<TubePointType>::iterator pnt; 
  typename ImageType::SizeType size = m_DataMask->GetLargestPossibleRegion().GetSize();
  typename ImageType::IndexType index;

  for(pnt = tube->GetPoints().begin(); pnt != tube->GetPoints().end(); pnt++)
    {
    x = (int)((*pnt).GetPosition()[0]+0.5);
    y = (int)((*pnt).GetPosition()[1]+0.5);
    z = (int)((*pnt).GetPosition()[2]+0.5);

    index[0]=x;
    index[1]=y;
    index[2]=z;

    if(x>=0 && x< size[0] &&
      y>=0 && y< size[1] &&
      z>=0 && z<size[2])
      {
      m_DataMask->SetPixel(index, (PixelType)(m_TubeID+(m_TubePointCount/10000.0)));
      r = (*pnt).GetRadius();
      r *= 1.25;
      if(r>1 &&
        x-r>=0 && x+r<size[0] &&
        y-r>=0 && y+r<size[1]  &&
        z-r>=0 && z+r<size[2] )
        {
        for(i=-(int)r; i<=(int)r; i++)
          {
          for(j=-(int)r; j<=(int)r; j++)
            {
            for(k=-(int)r; k<=(int)r; k++)
              {
              index[0]+=k;
              index[0]+=j;
              index[0]+=i;
              m_DataMask->SetPixel(index, (PixelType)(m_TubeID+(m_TubePointCount/10000.0)));
              }
            }
          }
        }
      }  
    m_TubePointCount++;
    }
  return true; 
}

/**
 * Draw tube */
template<class TInputImage>
template<class TDrawMask>
void
RidgeExtractor<TInputImage>
::DrawTube(TDrawMask * drawMask, TubeType * tube)
{
  typedef typename TDrawMask::PixelType DrawPixelType;

  m_TubeID = tube->GetId();
  m_TubePointCount = 0;
  int x, y, z, i, j, k;
  double r;
  std::vector<TubePointType>::iterator pnt; 
  typename TDrawMask::SizeType size = drawMask->GetLargestPossibleRegion().GetSize();
  typename TDrawMask::IndexType index;

  for(pnt = tube->GetPoints().begin(); pnt != tube->GetPoints().end(); pnt++)
    {
    x = (int)((*pnt).GetPosition()[0]+0.5);
    y = (int)((*pnt).GetPosition()[1]+0.5);
    z = (int)((*pnt).GetPosition()[2]+0.5);

    index[0]=x;
    index[1]=y;
    index[2]=z;

    if(x>=0 && x< size[0] &&
      y>=0 && y< size[1] &&
      z>=0 && z<size[2])
      {
      m_DataMask->SetPixel(index, (DrawPixelType)(m_TubeID+(m_TubePointCount/10000.0)));
      r = (*pnt).GetRadius();
      if(r>1 &&
        x-r>=0 && x+r<size[0] &&
        y-r>=0 && y+r<size[1]  &&
        z-r>=0 && z+r<size[2] )
        {
        for(i=-(int)r; i<=(int)r; i++)
          {
          for(j=-(int)r; j<=(int)r; j++)
            {
            for(k=-(int)r; k<=(int)r; k++)
              {
              index[0]+=k;
              index[0]+=j;
              index[0]+=i;
              m_DataMask->SetPixel(index, (DrawPixelType)(m_TubeID+(m_TubePointCount/10000.0)));
              }
            }
          }
        }
      }  
    m_TubePointCount++;
    }
}

/** Set the idle call back */
template<class TInputImage>
void
RidgeExtractor<TInputImage>
::IdleCallBack(bool (*idleCallBack)())
{
  m_IdleCallBack = idleCallBack;
}

/** Set the status callback  */
template<class TInputImage>
void
RidgeExtractor<TInputImage>
::StatusCallBack(void (*statusCallBack)(char *, char *, int))
{
  m_StatusCallBack = statusCallBack;
}


}; // end namespace itk

#endif
