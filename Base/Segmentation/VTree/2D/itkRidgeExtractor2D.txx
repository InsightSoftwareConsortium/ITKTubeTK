/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkRidgeExtractor2D.txx,v $
  Language:  C++
  Date:      $Date: 2006/05/19 11:58:09 $
  Version:   $Revision: 1.38 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRidgeExtractor2D_txx
#define __itkRidgeExtractor2D_txx

#include "itkRidgeExtractor2D.h"
#include "itkMatrixMath.h"
#include <itkImageRegionIterator.h>
#include <itkMinimumMaximumImageFilter.h>

namespace itk
{
 
template <class TInputImage>
class RidgeExtractor2DSplineValue : public UserFunc<vnl_vector<int> *, double>
{
public:
  RidgeExtractor2DSplineValue(RidgeExtractor2D<TInputImage> * newRidgeExtractor2D)
  {
    m_Ridge = newRidgeExtractor2D;
  };
  double value(vnl_vector<int>* x)
  {
    return m_Ridge->Intensity(x);
  };
protected:
  RidgeExtractor2D<TInputImage> * m_Ridge;
};


/**
 * Constructor */
template<class TInputImage>
RidgeExtractor2D<TInputImage>
::RidgeExtractor2D()
{
 
  m_Debug = false;
  VERBOSE = true;
   
  m_DataOp = BlurImageFunction<ImageType>::New();
  m_DataOp->SetScale(2); // 1.5
  m_DataOp->SetExtent(3.1); // 3

  SetModeMR(true);
  
  m_XD = new VectorType(ImageDimension);

  m_AutoScale = true;
  m_AutoScaleUsed = 2; // 1.5
  m_DynamicScale = true;
  m_RadiusExtractor = NULL;

  m_RecoveryMax = 4; 
  m_SplineValueFunc = new RidgeExtractor2DSplineValue<TInputImage>(this);
  m_DataSpline = new SplineND(ImageDimension, m_SplineValueFunc, &m_DataSpline1D, &m_DataSplineOpt);
  
  m_DataSplineOpt.searchForMin(true); //bug !!! to see !!!
  
  m_XH= new VnlMatrixType(ImageDimension,ImageDimension);
  m_XHEVect= new VnlMatrixType(ImageDimension,ImageDimension);

  this->m_IdleCallBack = NULL;
  this->m_StatusCallBack = NULL;
   
  m_Tube = NULL; 
 
}

/**
 * Destructor */
template<class TInputImage>
RidgeExtractor2D<TInputImage>
::~RidgeExtractor2D()
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
  
  if(m_XD != NULL)
  {    
    delete m_XD;
  }
  m_DataSpline = NULL;

}

/**
 * Set the input image */
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::SetInputImage(ImagePointer inputImage )
{
  m_Image = TInputImage::New();
  m_Image = inputImage;
  
  if(m_Image.IsNotNull()) 
  {
  
    m_DataOp->SetInputImage(inputImage);
   
    
    if(m_Debug) 
    {
     
      Point<float,ImageDimension> x;
      x[0] = 255;
      x[1] = 255;

      vnl_vector<double> vectX(2);
      vectX(0)=x[0];
      vectX(1)=x[1];
      
      VnlIntVectorType ivectX(2);
      ivectX[0]=(int) x[0];
      ivectX[1]=(int) x[1];

      std::cout << "Ridge: SetInputImage: m_DataOp->Evaluate[255,255] = " 
               << m_DataOp->Evaluate(x) << std::endl;
      std::cout << "Ridge: SetInputImage: m_DataOp->EvaluateIsoI[255,255] = " 
                 << m_DataOp->EvaluateIsoI(x) << std::endl;
      std::cout << "Ridge: SetInputImage: intensity at [255,255] = " 
               << this->Intensity(&ivectX) << std::endl;
      std::cout << "Ridge: SetInputImage: m_DataSpline[255,255] = " 
                << m_DataSpline->value(vectX) << std::endl;
    }

      
    typedef MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
    typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
    minMaxFilter->SetInput(m_Image);
    minMaxFilter->Update();
    m_DataMin = minMaxFilter->GetMinimum();
    m_DataMax = minMaxFilter->GetMaximum();
    m_DataRange = m_DataMax-m_DataMin;

    if(m_Debug)
    {
      std::cout << "Minimum = " << m_DataMin << std::endl;
      std::cout << "Maximum = " << m_DataMax << std::endl;
      std::cout << "Data Range = " << m_DataRange << std::endl;
    }

    /** Test */
    Point<float,TInputImage::ImageDimension> point;
    point[0]=50;
    point[1]=50;
    
   // double tf = m_DataOp->EvaluateIsoI(point); // ValueISOI
 
    unsigned int size[ImageDimension];
    for(unsigned int i=0;i<ImageDimension;i++)
    {
      size[i]= m_Image->GetLargestPossibleRegion().GetSize()[i];
      m_ExtractBoundMax[i] = size[i]-1;
      m_ExtractBoundMin[i] = 0;
    }

    m_DataSpline->xMin(m_ExtractBoundMin.Get_vnl_vector());
    m_DataSpline->xMax(m_ExtractBoundMax.Get_vnl_vector());

    /** Allocate the mask image */
    m_DataMask = ImageType::New();

    RegionType region;
    SizeType imSize;
    for(unsigned int i=0;i<ImageDimension;i++)
    {
      imSize[i]=size[i];
    }
    region.SetSize(imSize);
    region.SetIndex(m_Image->GetLargestPossibleRegion().GetIndex());
    m_DataMask->SetLargestPossibleRegion( region );
    m_DataMask->SetBufferedRegion( region );
    m_DataMask->SetRequestedRegion( region );
  
    m_DataMask->SetOrigin(m_Image->GetOrigin());
    m_DataMask->SetSpacing(m_Image->GetSpacing());
  
    m_DataMask->Allocate();

    /** Fill the mask image with zeros */
    itk::ImageRegionIterator<ImageType> 
      it(m_DataMask,m_DataMask->GetLargestPossibleRegion());

    it.GoToBegin();
    while(!it.IsAtEnd())
    {
      it.Set(0);
      ++it;
    }

  } // end Image == NULL

}


/**
 * Set Mode CT */
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::SetModeCT(bool modeCT)
{
  m_ModeMR = !modeCT;
  SetModeMR(m_ModeMR);
}

 
/**
 * Set Mode MR */
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::SetModeMR(bool modeMR)
{
  m_ModeMR = modeMR;
  if(modeMR)
  {
    m_ModeMR = true;
    m_StepX = 0.2;
    m_ThreshT = 0.75;
    m_ThreshX = 1.0;
    m_ThreshEV = -0.002;
    m_ThreshP2Q2 = 0.001;        // RW 0.005 // TMI 0.005 // near zero = harder
    m_ThreshP2Q2Start = 0.01;    // RW 0.01  //
    m_ThreshEVRatio = 0.2;      // RW 0.18  // TMI 0.2 // near 1 = harder
    m_ThreshEVRatioStart = 0.1;
    m_DataSplineOpt.tolerance(0.001);        // TMI 0.001
    m_DataSplineOpt.xStep(0.15);            // TMI 0.15
  }
  else
  {
    m_ModeMR = false;
    m_StepX = 0.2;
    m_ThreshT = 0.75;
    m_ThreshX = 1.0;
    m_ThreshP2Q2 = 0.04;           // near zero = harder
    m_ThreshP2Q2Start = 0.08;
    m_ThreshEV = -0.004;
    m_ThreshEVRatio = -0.225;          // near -1 = harder
    m_ThreshEVRatioStart = -0.2;
  }
}


/**
 * Set Mode Retina */
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::SetModeRetina(bool itkNotUsed( modeRetina ))
{
  m_StepX = 0.2;
  m_ThreshT = 0.75;
  m_ThreshX = 1.0;
  m_ThreshEV = -0.002;
  m_ThreshP2Q2 = 0.04;        // near zero = harder
  m_ThreshP2Q2Start = 0.08;    
  m_ThreshEVRatio = -0.225;       // near 1 = harder
  m_ThreshEVRatioStart = -0.2;
  m_DataSplineOpt.tolerance(0.001);       
  m_DataSplineOpt.xStep(0.15);        
 
  m_DataOp->SetScale(2);//0.5//2.5, 4 for Diaphragm, 2 for retina 
  m_DataOp->SetExtent(3);//3
   
  m_AutoScale = true;
  m_AutoScaleUsed = true; 
  m_DynamicScale = true;
  m_DynamicScaleUsed = 0.5;
}


/**
 * Compute Tangents */
template<class TInputImage>
bool
RidgeExtractor2D<TInputImage>
::m_CalcTangents(TubeType * newTube)
{
  if(newTube->GetPoints().size() == 0)
  {
    return true;
  }   
 
  VectorType tangent;
  tangent.Fill(0.0);

  if(newTube->GetPoints().size() == 1) 
  {
   (*newTube->GetPoints().begin()).SetTangent(tangent);
    return true;
  }
   
  std::vector<TubePointType>::iterator i, j, k;
  k = newTube->GetPoints().end();
  k--;
   
  for(i=newTube->GetPoints().begin(); i!=k; i++) 
  {
    j = i;
    j++;
    tangent = (*j).GetPosition() - (*i).GetPosition();
    (*i).SetTangent(tangent);
  }
   
  (*k).SetTangent(tangent);
   
  return true;
}

/**
 * Set Data Min value */
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::SetDataMin(double dataMin)
{
  m_DataMin = dataMin;
  m_DataRange = m_DataMax-m_DataMin;
}

/**
 * Set Data Min value */
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::SetDataMax(double dataMax)
{
  m_DataMax = dataMax;
  m_DataRange = m_DataMax-m_DataMin;
}

/**
 * Set Exctract bound min*/
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::SetExtractBoundMin(IntVectorType* extractBoundMin)
{
  m_ExtractBoundMin = * extractBoundMin;
}

/**
 * Get Exctract bound min*/
template<class TInputImage>
typename RidgeExtractor2D<TInputImage>::IntVectorType*
RidgeExtractor2D<TInputImage>
::GetExtractBoundMin(void)
{
  return & m_ExtractBoundMin;
}
  
/**
 * Set Exctract bound max*/
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::SetExtractBoundMax(IntVectorType* extractBoundMax)
{
  m_ExtractBoundMax = * extractBoundMax;
}

/**
 * Get Exctract bound mix*/
template<class TInputImage>
typename RidgeExtractor2D<TInputImage>::IntVectorType*
RidgeExtractor2D<TInputImage>
::GetExtractBoundMax(void)
{
  return & m_ExtractBoundMax;
}


/**
 * Set the scale */
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::SetScale(double scale)
{
  m_DataSpline->newData(true);
  m_DataOp->SetScale(scale);
}

/**
 * Get the scale */
template<class TInputImage>
double
RidgeExtractor2D<TInputImage>
::GetScale(void)
{
  return m_DataOp->GetScale();
}

/**
 * Set the extent */
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::SetExtent(double extent)
{
  m_DataSpline->newData(true);
  m_DataOp->SetExtent(extent);
}


/**
 * Get the extent */
template<class TInputImage>
double
RidgeExtractor2D<TInputImage>
::GetExtent(void)
{
  return m_DataOp->GetExtent();
}

/**
 * Get the data spline */
template<class TInputImage>
SplineND* 
RidgeExtractor2D<TInputImage>   
::GetDataSpline(void)
{
  return m_DataSpline;
}
  
/**
 * Get the data spline 1D*/
template<class TInputImage>
Spline1D* 
RidgeExtractor2D<TInputImage>
::GetDataSpline1D(void)
{
  return & m_DataSpline1D;
}
  
/**
 * Get the data spline optimizer*/
template<class TInputImage>
Optimizer1D* 
RidgeExtractor2D<TInputImage>
::GetDataSplineOptimizer(void)
{
  return & m_DataSplineOpt;
}


/**
 * Set the autoscale */
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::SetAutoScale(bool autoScale)
{
  if(m_RadiusExtractor)
  {
    m_AutoScale = autoScale;
  }
}


/**
 * Set the dynamic scale */
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::SetDynamicScale(bool dynamicScale)
{
  if(m_RadiusExtractor)
  {
    this->m_DynamicAutoScale = dynamicScale;
  }
}


/**
 * Set the radius extractor */
template<class TInputImage>
void 
RidgeExtractor2D<TInputImage>
::SetRadiusExtractor(RadiusExtractor2D<TInputImage> * radiusExtractor)
{
  m_RadiusExtractor = radiusExtractor;
}


/**
 * Return the intensity */
template<class TInputImage>
double 
RidgeExtractor2D<TInputImage>
::Intensity(VnlIntVectorType * x)
{
  Point<float,TInputImage::ImageDimension> point;
  for(unsigned int i=0;i<ImageDimension;i++)
  {
    point[i]=(*x)(i);
  }
  
  double tf;
  if(GetScale() > 0)
  {
    tf = (m_DataOp->EvaluateIsoI(point)-m_DataMin)/m_DataRange; // ValueISOI
  }
  else
  {
    Index<2> index;
    index[0]=static_cast<long int>(point[0]);
    index[1]=static_cast<long int>(point[1]);
    tf = (m_Image->GetPixel(index)-m_DataMin)/m_DataRange;

  }
  
  if(!m_ExtractRidge)
  {
    tf = 1-tf;
  }
  if(tf<0)
  {
    tf = 0;
  }
  else if(tf>1)
  {
    tf = 1;
  }
  
  if(m_Debug)
     std::cout << "   Ridge: Intensity: " << (*x)(0) << ", " << (*x)(1)<< " = " << tf << std::endl;
   
  return tf;
}

/**
 * Ridgeness
 */
template<class TInputImage>
double 
RidgeExtractor2D<TInputImage>
::Ridgeness(VectorType * x)
{    
 
  m_X = x;
   
  m_XVal = m_DataSpline->valueJet((*m_X).Get_vnl_vector(),(*m_XD).Get_vnl_vector(), *m_XH);
  
  unsigned int flat = 0;
  for(unsigned int i=0;i<ImageDimension;i++)
    {
    if(fabs((*m_XD)[i]) < 1e-10) 
      {
      flat++;
      }
    }

  if(flat == ImageDimension)
    {
    std::cout << "Ridgeness: Warning: Encountered flat region" << std::endl;

//    m_XH = 0;
 //  m_XHEVect = 0;
    m_XHEVal.Fill(0);
 //   m_XP = m_XQ = m_XR = 1;
    return 2;
    }

  Eigen(*m_XH, *m_XHEVect, m_XHEVal.Get_vnl_vector(),false);

  m_XP = (*m_XHEVect).get(0,0)*(*m_XD)[0] + (*m_XHEVect).get(1,0)*(*m_XD)[1];

  m_XD->Normalize();

  if(m_Debug) 
    {
    std::cout << "m_X " 
              << (*m_X)[0] << ", " << (*m_X)[1]  << std::endl;
    std::cout << "m_XVal " << m_XVal << std::endl;
    std::cout << "m_XP " << m_XP << std::endl;
    std::cout << "m_XQ " << m_XQ << std::endl;
    std::cout << "m_XR " << m_XR << std::endl;
    std::cout << "m_XD " 
              << (*m_XD)[0] << ", " << (*m_XD)[1] << std::endl;
    std::cout << "m_XHEVal " 
              << (m_XHEVal)[0] << ", " << (m_XHEVal)[1] << std::endl;
    std::cout << "m_XHEVect(0) " 
              << (*m_XHEVect)(0,0) << ", " << (*m_XHEVect)(1,0) << ", " 
              << std::endl;
    std::cout << "m_XHEVect(1) " << (*m_XHEVect)(0,1) << ", " 
              << (*m_XHEVect)(1,1) << ", "  << std::endl;
    std::cout << "m_XHEVect(0)*m_XHEVect(1) = " 
              << (*m_XHEVect)(0,0)*(*m_XHEVect)(0,1) +
                 (*m_XHEVect)(1,0)*(*m_XHEVect)(1,1)
              << std::endl;
    } 

  return m_XP*m_XP;

}

/**
 * Traverse one way
 */
template<class TInputImage>
typename RidgeExtractor2D<TInputImage>::TubePointer 
RidgeExtractor2D<TInputImage>::TraverseOneWay(VectorType * newX,VectorType * newT,VectorType  * newN,int dir)
// template<class TInputImage> void RidgeExtractor2D<TInputImage>::TraverseOneWay(VectorType * newX,VectorType * newT,VectorType  * newN,int dir)

{

  double lVal;
  VectorType lX;
  VectorType lT;
  VectorType lN;
  VectorType lNTEVal;
  VectorType lStepDir;
  VectorType lSearchDir;
  
  VectorType pX;
  VectorType pT;
  VectorType pN;
  VectorType pStepDir;
  VectorType pSearchDir;


  lX = (*newX);
  lT = (*newT);
  lN = (*newN);
  lStepDir = (*newT);
  lSearchDir = (*newN);
  
  pX = lX;
  pT = lT;
  pN = lN;
  pStepDir = lStepDir;
  pSearchDir = lSearchDir;
  
  VectorType prod(ImageDimension);
  VectorType tV(ImageDimension);
  
  Index<2> index;
  index[0] = static_cast<long int>( (*newX)[0]+0.5 );
  index[1] = static_cast<long int>( (*newX)[1]+0.5 );
  PixelType value = m_DataMask->GetPixel(index);

  if(value != 0 &&
     value != m_TubeID)
  {
    std::cout << "m_DataMask->GetPixel(index) != 0" << std::endl;
    return m_Tube;
  }
  else
  {
    m_DataMask->SetPixel(index,m_TubeID+(m_TubePointCount/10000.0));
  }

  double p2q2;
  double iScale0 = GetScale();
   
  int pSize = m_Tube->GetPoints().size();
   
  double stepFactor0 = 1;
  double stepFactor = stepFactor0;
   
  int i, j, k;
  double tf;
  int recovery = 0;
  int prevRecoveryPoint = -10;
  while(recovery < 4) 
  {
    if(recovery > 0) 
    {
      std::cout << "Attempting recovery : " << recovery << std::endl;            
      switch(recovery)
      {
        default:
   
        case 1:
               if(prevRecoveryPoint+10 <= m_TubePointCount)
               {
                 recovery++;
               }
               stepFactor = 0.8 * stepFactor0;
               //cDataSplineOpt.tolerance(0.000001);
               SetScale(iScale0 * 1.1);
               break;
       case 2:
               stepFactor = 1.5 * stepFactor0;
               SetScale(iScale0 * 1.2);
               break;
       case 3:
               stepFactor = 0.5 * stepFactor0;
               SetScale(iScale0 * 0.9);
               break;
      }

      std::cout << "Point = " << m_TubePointCount << ": Recovery: new scale = " << GetScale() << std::endl;
      prevRecoveryPoint = m_TubePointCount;
                
      lX = pX;
      lT = pT;
      lN = pN;
      lStepDir = pStepDir;
      lSearchDir = pSearchDir;
    }
    else 
    {
      if(GetScale() != iScale0)
      {
        SetScale(iScale0);
      } 
         //cDataSplineOpt.tolerance(0.0001);
                                  
      if(m_ThreshT>0 && fabs(dot_product(lStepDir.Get_vnl_vector(), pStepDir.Get_vnl_vector()))<1-0.5*(1-m_ThreshT)) //vDotProd works in 2D too et fait la somme de la multipication termet a terme
      {  
        std::cout << "Curving" << std::endl;
        stepFactor *= (float)0.75;
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
      
    lX.Set_vnl_vector(ComputeLineStep(lX.Get_vnl_vector(), m_StepX*stepFactor, lStepDir.Get_vnl_vector()));


    if(m_Debug) 
      {
      std::cout << "Ridge: TraverseOW: lStepDir = " << lStepDir[0] << ", " << lStepDir[1] << std::endl;
      std::cout << "Ridge: TraverseOW: lX0 = " << lX[0] << ", " << lX[1] << std::endl;
      std::cout << "Ridge: TraverseOW: lSearchDir1 = " << lSearchDir[0] << ", " << lSearchDir[1] << std::endl;
      }

    SplineND::VectorRefType lXVnlRef = lX.Get_vnl_vector();
    SplineND::VectorType LSearchDirVnl = lSearchDir.Get_vnl_vector();
    
  //  if(!m_DataSpline->extreme(lX.Get_vnl_vector(), &lVal, lSearchDir.Get_vnl_vector()))
    if(!m_DataSpline->extreme(lXVnlRef, &lVal, LSearchDirVnl))
      {
      std::cout << "*** Ridge terminated: Local max not found" << std::endl;
      recovery++;
      continue;
      }
    
    if(m_Debug)
      {
      std::cout << "Ridge: TraverseOW: lXExtreme = " << lX[0] << ", " << lX[1] << std::endl;
      }
      
 
    p2q2 = Ridgeness(&lX);  

    if ( p2q2 == 2 ) 
      {
      std::cout << "*** Ridge terminated: flatness" << std::endl;
      recovery += 4 ;
      continue;
      }
   // k = 1;
      
    VectorType prod;    
    VectorType tV;
      
    for(i=0; i<=1; i++) 
      {
      for(j=0; j<=1; j++) 
        {
        tV[j] = (*m_XHEVect).get(j, i);
        }
      prod[i] = fabs(dot_product(lStepDir.Get_vnl_vector(), tV.Get_vnl_vector()));
      }
 
    tf = prod[0];
    
    k=2;
    if(prod[1]>tf)
    {
      k = 2;
      tf = prod[1];
    }
        
    if(k != 2) 
    {
      std::cout << "Mixing things up: Chosen t=evect #" <<k<< " tf=" << tf <<std::endl;
    }   
     
    j = (k)%2;
   
    for(i=0; i<=1; i++)
    {
      lN[i] = (*m_XHEVect).get(i, j);
    }    
    lNTEVal[1] = (m_XHEVal)[j];
    lSearchDir = lN;
         
    for(i=0; i<=1; i++)
      {
      lT[i] = (*m_XHEVect).get(i, k-1);
      }
         
    tf = dot_product(pT.Get_vnl_vector(), lT.Get_vnl_vector());     
    if(tf<0) 
      {
      lT[0] = lT[0] * -1;
      lT[1] = lT[1] * -1;
      }
         
    lStepDir[0] = (pStepDir[0] * 0.25) + (lT[0] * 0.75);
    lStepDir[1] = (pStepDir[1] * 0.25) + (lT[1] * 0.75);
       
    tf = dot_product(lStepDir.Get_vnl_vector(), pStepDir.Get_vnl_vector());

    if(m_ThreshT>0 && fabs(tf)<m_ThreshT)
    {
      std::cout << "*** Ridge terminated: Rapid change in step direction " 
                << "(" << fabs(tf) << ")" << std::endl;
      recovery++;     
      continue;
    }
         
    tf = sqrt(ComputeEuclideanDistanceVector(lX.Get_vnl_vector(), pX.Get_vnl_vector()));
    if(m_ThreshX>0 && tf>m_ThreshX) 
    {
      std::cout << "*** Ridge terminated: Rapid change in spatial location "
                << "(" << tf << ")" << std::endl;
      recovery++;
      continue;
    }
   
    if(p2q2>m_ThreshP2Q2) 
    {
      std::cout << "*** Ridge terminated: Local max not a ridge point " << "(p2q2 = " << p2q2 << ")" << std::endl;
         
      tf = 1/(lNTEVal[0]);
      std::cout << "       EV Ratio = " << tf << std::endl;
      recovery++;
      continue;
    }
         
    
    if((int)(lX[0]+0.5)<(m_ExtractBoundMin)[0] || (int)(lX[1]+0.5)>(m_ExtractBoundMax)[0] ||
       (int)(lX[1]+0.5)<(m_ExtractBoundMin)[1] || (int)(lX[2]+0.5)>(m_ExtractBoundMax)[1]) 
    {
      std::cout << "*** Ridge terminated: Exited extraction bounds" << std::endl;
      break;
    }
         
    index[0]= static_cast <unsigned int> ( lX[0]+0.5 );
    index[1]= static_cast <unsigned int> ( lX[1]+0.5 );

    tf = m_DataMask->GetPixel(index);
    
    if(tf != 0) 
    {
      if((int)tf != m_TubeID || m_TubePointCount-((tf-(int)tf)*10000)>5/m_StepX) 
      {
        std::cout << "*** Ridge terminated: Revisited voxel" << std::endl;
        break;
      }
    }
    else
    {
      m_DataMask->SetPixel(index,m_TubeID+(m_TubePointCount/10000.0));
    }   
         
    if(m_TubePointCount%50==0)
    {
      if(this->m_StatusCallBack) 
      {
        char st[80];
        sprintf(st, "Point #%d", m_TubePointCount);
        this->m_StatusCallBack(NULL, st, 0);
      }
      else
      {
//        std::cout << "Adding point " << m_TubePointCount << " = " << lX(1) << ", " << lX(2) << std::endl;
      }
    }
     
    TubePointType pnt;
    pnt.SetID(m_TubePointCount);
    pnt.SetPosition(lX[0],lX[1]);
    pnt.SetRidgeness(p2q2);
    //pnt->SetColor(&cTubeColor);

    if(dir == 1)
      {
      //m_Tube->GetPoints().push_front(pnt);
      m_TubePointList.push_front(pnt);
      }
    else
      {
      //m_Tube->GetPoints().push_back(pnt);
      m_TubePointList.push_back(pnt);
      }   
    
    m_TubePointCount++;
         
    recovery = 0;
         
    if(m_TubePointCount/25.0 == m_TubePointCount/25) 
    {
      if(this->m_IdleCallBack!=NULL)
      {
        this->m_IdleCallBack();
      }
      if(m_DynamicScale && m_TubePointCount/50.0 == m_TubePointCount/50)
      {
        if(m_RadiusExtractor != NULL) 
        {       
          TubePointType tmpPoint;
          tmpPoint.SetPosition(lX[0],lX[1]);
          tmpPoint.SetTangent(lT[0],lT[1]);
          tmpPoint.SetRadius(GetScale());
          

          m_RadiusExtractor->CalcOptimalScale(tmpPoint);
          m_DynamicScaleUsed = (m_DynamicScaleUsed+tmpPoint.GetRadius()*0.66)/2.0;
          SetScale(m_DynamicScaleUsed);        
          if(m_DynamicScaleUsed<1.0)
          {
            m_DynamicScaleUsed = 1.0;
          }   
          if(this->m_StatusCallBack) 
          {
            char s[80];
            sprintf(s, "Extract: Ridge: DS = %1.1f", m_DynamicScaleUsed);
            this->m_StatusCallBack(s, NULL, 0);
          }
        }
        else
        {
          std::cout << "Dynamic Scale = " << m_DynamicScaleUsed << std::endl;
          SetScale(m_DynamicScaleUsed);
        }
      }
    }  
  }
    
  std::cout << "*** Ridge terminated: Cannot recover" << std::endl;
  std::cout << "    Length = " << m_Tube->GetPoints().size()-pSize << std::endl;

  SetScale(iScale0);

  return m_Tube;
}
  

/** Compute the local ridge */
template<class TInputImage>
bool
RidgeExtractor2D<TInputImage>
::LocalRidge(VectorType * newX)
{
  Index<2> index;
  index[0] = static_cast<long int>( (*newX)[0]+0.5 );
  index[1] = static_cast<long int>( (*newX)[1]+0.5 );
 
  if(m_DataMask->GetPixel(index) != 0)
  {      
    if(this->m_StatusCallBack) 
    {
      this->m_StatusCallBack(NULL, "Revisited voxel", 0);
      return false;
    }
    else 
    {
      std::cout << "Revisited voxel" << std::endl;
      return false;
    }
  }
  
  if(this->m_StatusCallBack)
  {
    this->m_StatusCallBack(NULL, "SEARCHING", 0);
  }

  int i;  
  double val;

  VnlVectorType lN(ImageDimension);
  //VnlMatrixType lN(ImageDimension,ImageDimension-1);
  Ridgeness(newX);

  
  for(i=0; i<=1; i++) 
  {
    lN[i] = (*m_XD)[i];
  }
  
  m_DataSpline->extreme((*newX).Get_vnl_vector(), &val, lN);
  
  double p2q2;
  for(i=0; i<=1; i++) 
  {
    lN[i] = (*m_XHEVect).get(i, 0);
  }

  m_DataSpline->extreme((*newX).Get_vnl_vector(), &val, lN);
   
 
  for(i=0; i<=1; i++) 
  {
    lN[i] = (*m_XHEVect).get(i, 0);
  }
   
  m_DataSpline->extreme((*newX).Get_vnl_vector(), &val, lN);
   

  index[0] =  static_cast<long int>( (*newX)[0]+0.5 );
  index[1] =  static_cast<long int>( (*newX)[1]+0.5 );


  if(m_DataMask->GetPixel(index) != 0)
  {
    if(this->m_StatusCallBack) 
    {
      this->m_StatusCallBack(NULL, "Revisited voxel", 0);
      return false;
    }
    else 
    {
      std::cout << "Revisited voxel" << std::endl;
      return false;
    }
  }
  

  for(i=0; i<=1; i++) 
  {
    lN(i) = (*m_XHEVect).get(i, 0);
  }

  if(!m_DataSpline->extreme((*newX).Get_vnl_vector(), &val, lN)) 
  {
    if(this->m_StatusCallBack)
    {
      this->m_StatusCallBack(NULL, "Not Found: No local max", 0);
    }
    std::cout << "*** Ridge not found: Local max not found" << std::endl;
    return false;
  }

  index[0] =  static_cast<long int>( (*newX)[0]+0.5 );
  index[1] = static_cast<long int>( (*newX)[1]+0.5 );
  if(m_DataMask->GetPixel(index) != 0)
  {
    if(this->m_StatusCallBack) 
    {
      this->m_StatusCallBack(NULL, "Revisited voxel", 0);
      return false;
    }
    else 
    {
      std::cout << "Revisited voxel" << std::endl;
      return false;
    }
  }
   
  p2q2 = Ridgeness(newX);

  if(p2q2 > m_ThreshP2Q2Start) 
  {
    if(this->m_StatusCallBack)
    {
      this->m_StatusCallBack(NULL, "Not Found: P2+Q2 != 0", 0);
    }
    std::cout << "*** Ridge not found: P2+Q2 (" << p2q2 << ") not near zero" << std::endl;
    std::cout << "       EV Ratio = " << (m_XHEVal)[1]/(float)sqrt(((m_XHEVal)[0]*(m_XHEVal)[0]+(m_XHEVal)[1]*(m_XHEVal)[1])/2) << std::endl;
    return false;
  }
   
  return true;
}



/**
 * Extract a tube 
 */
template<class TInputImage>
typename RidgeExtractor2D<TInputImage>::TubePointer
RidgeExtractor2D<TInputImage>
::Extract(VectorType * newX, int tubeID)
{
  
  VectorType lX;
  lX = (*newX);
//  std::cout<<lX[0]<<std::endl;
  if(!LocalRidge(&lX))
  {
  
    if(m_Debug) std::cout << "LocalRidge fails !" << std::endl;
    return NULL;
  }
  std::cout << "*** Ridge found at " << lX[0] << ", " << lX[1]  << std::endl;

  m_TubePointList.clear();

  int i;
  VectorType lN;
  VectorType lT;
   
  double iScaleOld = GetScale();

  m_RadiusExtractor->SetRadius0(iScaleOld);
  
  if(m_AutoScale && m_RadiusExtractor != NULL) 
  {
    for(i=0; i<ImageDimension; i++)
    {
      lT[i] = (*m_XHEVect).get(i,1);
    } 
    TubePointType tmpPoint;
    tmpPoint.SetPosition(lX[0],lX[1]);
    tmpPoint.SetTangent(lT[0],lT[1]);
    tmpPoint.SetRadius(GetScale());
   
   
    if(!m_RadiusExtractor->CalcOptimalScale(tmpPoint,true)) 
    {
      if(m_Debug && this->m_StatusCallBack)
      {
        this->m_StatusCallBack("Extract: Ridge: AS = ?", "Error: Medial Max Not Found", 0);
      }
      m_AutoScaleUsed = GetScale();
    }
    else 
    {
      m_AutoScaleUsed = tmpPoint.GetRadius()*0.66;//(2*tmpPoint.GetR()+m_AutoScaleUsed)/3;
    }

    if(m_AutoScaleUsed<1.0)
    {
      m_AutoScaleUsed = 1.0;
    }
    m_DynamicScaleUsed = m_AutoScaleUsed;
    SetScale(m_AutoScaleUsed);

    //m_RadiusExtractor->SetRadius0(m_AutoScaleUsed);
   
    if(VERBOSE)
    {   
      std::cout << "AutoScale = " << GetScale() << std::endl;
    }

    lX = (*newX);
    if(!LocalRidge(&lX)) 
    {
      if(this->m_StatusCallBack)
      {
        this->m_StatusCallBack("AutoScale Failure", NULL, 0);
      }
      if(m_Debug)
      {
        std::cout << "RidgeExtractor2D:Extract(): AS Failure" << std::endl;
      }
      //m_AutoScaleUsed = iScaleOld;
      //SetScale(iScaleOld);
      return NULL;
    }
  }
  
  m_Tube = TubeType::New();
  m_TubeID = tubeID;
  m_Tube->SetId(tubeID);
  m_TubePointCount = 0;

  for(i=0; i<ImageDimension; i++) 
  {
    lN[i] = m_XHEVect->get(i,0);
    lT[i] = m_XHEVect->get(i,ImageDimension-1);
  }

  TraverseOneWay(&lX, &lT, &lN, 1);
  


  
  if(m_AutoScale)
  {
    SetScale(m_AutoScaleUsed);
  }
  else
  {
    SetScale(iScaleOld);
  }

  for(unsigned int i=0;i<ImageDimension;i++)
  {
    lT[i] = -1 * lT[i];
  }

  TraverseOneWay(&lX, &lT, &lN, -1);

   // Set the list of tubepoints
  std::list<TubePointType>::const_iterator ptIt = m_TubePointList.begin();
  while( ptIt != m_TubePointList.end())
    {
    m_Tube->GetPoints().push_back(*ptIt);
    ptIt++;
    }


  SetScale(iScaleOld);
     
  if(m_Tube->GetPoints().size()<10) 
  {
    if(this->m_StatusCallBack)
    {
      this->m_StatusCallBack("Extract: Ridge", "Too short", 0);
    }
    if(VERBOSE)
    {
      std::cout << "Extract: Ridge Too short" << std::endl;
    }
    DeleteTube(m_Tube);
    m_Tube = NULL;
    return m_Tube;
  }
  
  std::cout << "*** Extracted ridge of " << m_Tube->GetPoints().size() << " points." << std::endl;
  
  if(m_Tube.IsNotNull())
  {
    m_CalcTangents(m_Tube);
  } 
 
  if(this->m_StatusCallBack) 
  {
    char s[80];
    sprintf(s, "%d points", m_Tube->GetPoints().size());
    this->m_StatusCallBack("Extract: Ridge", s, 0);
  }   
  return m_Tube;
}


/**
 * Smooth a tube 
 */
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::SmoothTubeX(TubeType * tube, int h)
{
  int i, cnt;
  PointType avg;

  std::vector<TubePointType> & points = tube->GetPoints();
  std::vector<TubePointType>::iterator pnt, pntT;

  std::vector<TubePointType>::iterator begin = points.begin();
  std::vector<TubePointType>::iterator end = points.end();


  DeleteTube(tube);
  
  //  pnt = tube->GetPoints().begin();
     
  
  for(pnt = begin;
      pnt != end;
      pnt++)
  {
    
    cnt = 0;
    
    avg.Fill( 0);
   
    if(pnt != tube->GetPoints().begin())
    {
      pntT = pnt;
      for(i=0; i<h/2 && pntT!=tube->GetPoints().begin(); i++, pntT--)
      {
        avg[0] += (*pntT).GetPosition()[0];
        avg[1] += (*pntT).GetPosition()[1];
        cnt++;
      }
    }
    if(pnt != tube->GetPoints().end())
    {
      pntT = pnt;
      pntT++;
      for(i=0; i<h/2 && pntT!=tube->GetPoints().end(); i++, pntT++)
      {
        avg[0] += (*pntT).GetPosition()[0];
        avg[1] += (*pntT).GetPosition()[1];
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
 

/**
 * Delete a tube 
 */
template<class TInputImage>
bool
RidgeExtractor2D<TInputImage>
::DeleteTube(TubeType * tube)
{    
  int x, y;
  int j, k;
  double r;
  std::vector<TubePointType>::iterator pnt;

  for(pnt = tube->GetPoints().begin(); pnt != tube->GetPoints().end(); pnt++)
  {
    x = (int)((*pnt).GetPosition()[0]+0.5);
    y = (int)((*pnt).GetPosition()[1]+0.5);
    IndexType index;
    index[0]=x;
    index[1]=y; 
    m_DataMask->SetPixel(index,0);
    r = (*pnt).GetRadius();
    r /= 2.0;
    if(r>1 &&
       x-r>=0 && x+r<m_DataMask->GetLargestPossibleRegion().GetSize()[0] &&
       y-r>=0 && y+r<m_DataMask->GetLargestPossibleRegion().GetSize()[1]
       )
    {
      for(j=static_cast<int>(-r); j<=r; j++)
      {
        for(k=static_cast<int>(-r); k<=r; k++)
        {
          index[0]=x+k;
          index[1]=y+j;
          m_DataMask->SetPixel(index,0);
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
RidgeExtractor2D<TInputImage>
::AddTube(TubeType * tube)
{
  m_TubeID = tube->GetId();
  m_TubePointCount = 0;
  int x, y, j, k;
  double r;
  std::vector<TubePointType>::iterator pnt; 
  
  SizeType size = m_DataMask->GetLargestPossibleRegion().GetSize();
 
  IndexType index;

  for(pnt = tube->GetPoints().begin(); pnt != tube->GetPoints().end(); pnt++)
  {
    x = (int)((*pnt).GetPosition()[0]+0.5);
    y = (int)((*pnt).GetPosition()[1]+0.5);
    index[0]=x;
    index[1]=y;

      if(x>=0 && (unsigned long)x< size[0] &&
        y>=0 && (unsigned long)y< size[1])
      {
        m_DataMask->SetPixel(index,m_TubeID+(m_TubePointCount/10000.0));
        r = (*pnt).GetRadius();
        r /= 2.0;
        if(r>1 &&
           x-r>=0 && x+r<size[0] &&
           y-r>=0 && y+r<size[1]
          )
        {
          for(j=static_cast<int>(-r); j<=r; j++)
            for(k=static_cast<int>(-r); k<=r; k++)
            {
              index[0]+=k;
              index[1]+=j;
              m_DataMask->SetPixel(index,m_TubeID+(m_TubePointCount/10000.0));
            }
      }
    }

    m_TubePointCount++;
  }
  return true; 
}

/**
 * Set the idle call back */
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::IdleCallBack(bool (*idleCallBack)())
{
  this->m_IdleCallBack = idleCallBack;
}

/**
 * Set the status callback  */
template<class TInputImage>
void
RidgeExtractor2D<TInputImage>
::StatusCallBack(void (*statusCallBack)(char *, char *, int))
{
  this->m_StatusCallBack = statusCallBack;
}


}; // end namespace itk

#endif
