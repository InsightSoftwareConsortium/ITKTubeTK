/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkRadiusExtractor2D.txx,v $
  Language:  C++
  Date:      $Date: 2006/05/19 04:00:55 $
  Version:   $Revision: 1.22 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRadiusExtractor2D_txx
#define __itkRadiusExtractor2D_txx

#include "itkRadiusExtractor2D.h"
#include <itkMinimumMaximumImageFilter.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

namespace itk
{

/** Define the Medialness Function */
template <class T>
class RadiusExtractor2DMedialnessFunc : public UserFunc<double, double> 
{
   typedef itk::VesselTubeSpatialObject<2> TubeType;
   typedef typename TubeType::TubePointType TubePointType;
   std::list<TubePointType>* m_Tube;
   RadiusExtractor2D<T> * m_RadiusExtractor2D;

public:
   RadiusExtractor2DMedialnessFunc(RadiusExtractor2D<T> * newRadiusExtractor, std::list<TubePointType> * newTube)
   {
      m_Tube = newTube;
      m_RadiusExtractor2D = newRadiusExtractor;
   };
   double value(double x)
   {
     return m_RadiusExtractor2D->MedialnessAtKern(m_Tube, x);
   };

};

   
/**
 * Constructor */
template<class TInputImage>
RadiusExtractor2D<TInputImage>
::RadiusExtractor2D()
{
  m_Debug = false; 
  m_ModeMR = true;
   
  m_DataOp = BlurImageFunction<ImageType>::New();
  m_DataOp->SetScale(1.0);
  m_DataOp->SetExtent(1.1);

  m_NumRadiusPoints = 5;
  m_RadiusPointSpacing = 10;
  m_Radius0 = 1.0;
  m_ExtractRidge = true;
   
  m_ThreshWVal = 0.04;       // 0.015; larger = harder
  m_ThreshWValStart = 0.01;
   
  m_KernNumT = 0;
  double theta;
  for(theta=0; theta<PI-PI/8; theta+=(double)(PI/4)) 
  {
    m_KernSinT[m_KernNumT] = sin(theta);
    m_KernCosT[m_KernNumT] = cos(theta);
    m_KernNumT++;
  }
  m_KernN0 = new VectorType(ImageDimension);
  m_KernN1 = new VectorType(ImageDimension);
  for(unsigned int i=0;i<ImageDimension;i++)
  {
    (*m_KernN0)(i) = 0;
    (*m_KernN1)(i) = 0;
  }
  
  (*m_KernN0)(0) = 1;
  (*m_KernN1)(1) = 1;
  typedef std::list<TubePoint *>::iterator ListIteratorType;
 // m_IterPntArray = new typedef std::list<TubePoint *>::iterator[5000]; 
  m_IterPntArray = new ListIteratorType [5000];
  m_KernPntArray = new TubePoint *[5000]; 
  m_ArrayLen = 0;
  std::list<TubePoint *> InitList;
  
  for(int i=0; i<5000; i++)
    {
    m_IterPntArray[i] = InitList.end();
    m_KernPntArray[i] = NULL;
    }

  m_Kern = new std::list<TubePointType>;

  m_MedialnessAtKern = new RadiusExtractor2DMedialnessFunc<TInputImage>(this,m_Kern);
  m_MedialnessOpt.use(m_MedialnessAtKern);
  m_MedialnessOpt.xMin(0.4);
  m_MedialnessOpt.xMax(20.0);

  m_MedialnessOpt.tolerance(0.001); 
  m_MedialnessOpt.xStep(0.25);

  m_MedialnessOpt.searchForMin(false);
   
  m_IdleCallBack = NULL;
  m_StatusCallBack = NULL;
}

/** 
 * Destructor */
template<class TInputImage>
RadiusExtractor2D<TInputImage>
::~RadiusExtractor2D()
{
  if(m_MedialnessAtKern != NULL)
  {
    delete m_MedialnessAtKern;
  } 
 
  m_MedialnessAtKern = NULL;
   
  for(unsigned int i=0; i<(unsigned int)m_ArrayLen; i++)
  {
    if(m_KernPntArray[i] != NULL)
    {
      delete m_KernPntArray[i];
      m_KernPntArray[i] = NULL;
    }
  }
  
  m_ArrayLen = 0;
  delete [] m_IterPntArray;
  delete [] m_KernPntArray;
}


/** 
 * Set the scale factor */
template<class TInputImage>
void
RadiusExtractor2D<TInputImage>
::SetScale(double scale)
{
  m_Scale = scale;
  m_DataOp->SetScale(scale);
}


/**
 * Set the extent factor */
template<class TInputImage>
void
RadiusExtractor2D<TInputImage>
::SetExtent(double extent)
{
  m_Extent = extent;
  m_DataOp->SetExtent(extent);
}


/**
 * Set Mode CT */
template<class TInputImage>
void
RadiusExtractor2D<TInputImage>
::SetModeCT(bool modeCT)
{
  m_ModeMR = !modeCT;
}
 

/**
 * Set Mode Retina */
template<class TInputImage>
void
RadiusExtractor2D<TInputImage>
::SetModeRetina(bool itkNotUsed( modeRetina ))
{
  m_Debug = false; 
  m_ModeMR = true;
   
  m_DataOp->SetScale(0.3);
  m_DataOp->SetExtent(3);

  m_NumRadiusPoints = 5;
  m_RadiusPointSpacing = 3;
  m_Radius0 = 0.75;
  m_ExtractRidge = false;
   
  m_ThreshWVal = 0.01;       
  m_ThreshWValStart = 0.005;
 
  m_MedialnessOpt.xMin(0.4);
  m_MedialnessOpt.xMax(20.0);

  m_MedialnessOpt.tolerance(0.0001); 
  m_MedialnessOpt.xStep(0.05);

  m_MedialnessOpt.searchForMin(false);
}

/**
 * Set Radius Min */
template<class TInputImage>
void
RadiusExtractor2D<TInputImage>
::SetRadiusMin(double radiusMin)
{
  m_RadiusMin = radiusMin;
  m_MedialnessOpt.xMin(m_RadiusMin);
}


/**
 * Set Radius Max */
template<class TInputImage>
void
RadiusExtractor2D<TInputImage>
::SetRadiusMax(double radiusMax)
{
  m_RadiusMax = radiusMax;
  m_MedialnessOpt.xMax(m_RadiusMax);
}

/**
 * Get the medialness operator */
template<class TInputImage>
typename RadiusExtractor2D<TInputImage>::OptimizerType & 
RadiusExtractor2D<TInputImage>
::GetMedialnessOpt(void)
{
  return & m_MedialnessOpt;
}

/**
 * Set the input image */
template<class TInputImage>
void
RadiusExtractor2D<TInputImage>
::SetInputImage(ImagePointer inputImage )
{
  m_Image = inputImage;
  
  if(m_Image)
  {
    m_DataOp->SetInputImage(m_Image);
    typedef MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
    typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New();
    minMaxFilter->SetInput(m_Image);
    minMaxFilter->Update();
    m_DataMin = minMaxFilter->GetMinimum();
    m_DataMax = minMaxFilter->GetMaximum();

    if(m_Debug)
    {
      std::cout << "RadiusExtractor2D: SetInputImage: Minimum = " << m_DataMin << std::endl;
      std::cout << "RadiusExtractor2D: SetInputImage: Maximum = " << m_DataMax << std::endl;
    }
  }
}

/**
 * Compute the medialness at a given point */
template<class TInputImage>
double
RadiusExtractor2D<TInputImage>
::MedialnessAtPoint(TubePointType pnt,double pntR)
{
  VectorType y(2);
  VectorType n1(2);

  n1(0)=pnt.GetTangent()[1];
  n1(1)=-pnt.GetTangent()[0];
  n1.normalize();
 
  n1(0) /= m_Image->GetSpacing()[0];
  n1(1) /= m_Image->GetSpacing()[1];
    

  m_DataOp->SetExtent(2.1);
  m_DataOp->SetScale(2);
    
  int j;
  int c;
  double v, r;
  double kernPos[20], kernNeg[20];
  for(c=0; c<2; c++) 
  {
    if(c == 0)
    {
      r = 2*pntR/3;
    }
    else
    {
      r = 4*pntR/3;
    }
    j = 0;
       
    VectorType tp(ImageDimension);
    tp(0) = pnt.GetPosition()[0];
    tp(1) = pnt.GetPosition()[1];
    
    y=ComputeLineStep(tp, r, n1);
    
    if(m_Debug)
    {
      std::cout << "PointPosition = " << tp << std::endl;
      std::cout << "LineStep = " << y << std::endl;
    }

    if(y[0]>=0 && y[0]<m_Image->GetLargestPossibleRegion().GetSize()[0] &&
       y[1]>=0 && y[1]<m_Image->GetLargestPossibleRegion().GetSize()[1]
      ) 
    {
      Point<float,ImageDimension> point;
      
      for(unsigned int id=0;id<ImageDimension;id++)
      {
        point[id]=y(id);
      }
      v = (m_DataOp->Evaluate(point)-m_DataMin)/(m_DataMax-m_DataMin);
    if(m_Debug)
    {
      std::cout << "v = " << v << " m_DataMin = " << m_DataMin<<" m_DataMax = "<<m_DataMax << std::endl;
    }
      if(!m_ExtractRidge)
      { 
      v = 1-v;
      }
      if(v<0) 
      {
        v = 0;
      }
      if(v>1)
      {
        v = 1;
      }
      if(c == 0)
      { 
        kernPos[j] = v*v;
      }
      else
      {
        kernNeg[j] = v*v;
      }
    }
    j++;

    tp(0) = pnt.GetPosition()[0];
    tp(1) = pnt.GetPosition()[1];
    y= ComputeLineStep(tp, -r, n1);
           
    if(y[0]>=0 && y[0]<m_Image->GetLargestPossibleRegion().GetSize()[0] &&
       y[1]>=0 && y[1]<m_Image->GetLargestPossibleRegion().GetSize()[1] 
      ) 
    {   
      Point<float,ImageDimension> point;
      for(unsigned int id=0;id<ImageDimension;id++)
      {
        point[id]=y(id);
      }
      v = (m_DataOp->Evaluate(point)-m_DataMin)/(m_DataMax-m_DataMin);
      if(!m_ExtractRidge)
        {
        v = 1-v;
        }   
   
      if(v<0)
        {
        v = 0;
        }       
      if(v>1) 
        {
        v = 1;
        }
      if(c == 0)
        {
        kernPos[j] = v*v;
        }
      else
        {
        kernNeg[j] = v*v;
        }
    }
    j++;   
  }
    
  int k;
  double filtTot, kernTot[200];
 kernTot[0] = kernPos[0]-kernNeg[0];
  filtTot = kernTot[0];
  
   if(m_Debug)
      {
      for (int tempI =0 ; tempI < j ; tempI ++ )
        {
        std::cout<<tempI<<" kernPos "<< kernPos[tempI] <<" kernNeg " << kernNeg[tempI]<<std::endl;
        }
      }


  for(k=1; k<j; k++) 
  {
    kernTot[k] = kernPos[k]-kernNeg[k];
    filtTot += kernTot[k]; 
  }
 
 if ( j==0 )
    {
    filtTot = 0;
    }
  else
    {
    filtTot /= (j);
    }
  if (filtTot<-9999||filtTot>99e20)
    {
    filtTot = 0;
    }
  pnt.SetMedialness(filtTot); // Useless this line
  
 
  return filtTot;
}

/**
 * Compute the branchness at a given point
 */
template<class TInputImage>
double 
RadiusExtractor2D<TInputImage>
::BranchnessAtPoint(TubePointType pnt,double pntR)
{
  VectorType y(ImageDimension);
  VectorType n1(ImageDimension);

  n1(0)=(pnt.GetTangent())(1);
  n1(1)=-((pnt.GetTangent())(0));
  n1.normalize();
 
  n1(0) /= m_Image->GetSpacing()[0];
  n1(1) /= m_Image->GetSpacing()[1];  

  double kernPos;
  m_DataOp->SetScale(pntR);
  m_DataOp->SetExtent(4);
  y = pnt.GetPosition();
  
  Point<double,ImageDimension> point;   
  for(unsigned int id=0;id<ImageDimension;id++)
  {
    point[id]=y(id);
  }
  kernPos = (m_DataOp->Evaluate(point)-m_DataMin)/(m_DataMax-m_DataMin);

  if(!m_ExtractRidge)
  {
        kernPos = 1-kernPos;
  }
  
  kernPos = kernPos*kernPos;
  pntR *= 2.0;//normalement 2.0
  double kernNeg[20];
  double e = 6+(-1.5-log(pntR/4));
  if(e<2) {e = 2;}
  if(e>6) {e = 6;}
  m_DataOp->SetExtent(e);
  m_DataOp->SetScale(pntR/4);
  int j;
  double v;
  double theta;
  j = 0;

  for(theta=0; theta<PI-PI/16; theta+=(double)(PI/8)) 
  {
    y=ComputeLineStep(pnt.GetPosition(), pntR*(double)sin(theta), n1);
    y=ComputeLineStep(y, pntR*(double)cos(theta), n1);

    if(y[0]>=0 && y[0]<m_Image->GetLargestPossibleRegion().GetSize()[0] 
       && y[1]>=0 && y[1]<m_Image->GetLargestPossibleRegion().GetSize()[1]
      ) 
    { 
      Point<double,ImageDimension> point;   
      for(unsigned int id=0;id<ImageDimension;id++)
      {
        point[id]=y(id);
      }
      v = (m_DataOp->Evaluate(point)-m_DataMin)/(m_DataMax-m_DataMin);
      if(!m_ExtractRidge)
      {
        v = 1-v;
      }
      if(v<0) { v = 0; }
      if(v>1) { v = 1; }
      kernNeg[j] = v*v;
    }
    
    j++;

    y=ComputeLineStep(pnt.GetPosition(), -pntR*(double)sin(theta), n1);
    y=ComputeLineStep(y, -pntR*(double)cos(theta), n1);

    if(y[0]>=0 && y[0]<m_Image->GetLargestPossibleRegion().GetSize()[0] &&
       y[1]>=0 && y[1]<m_Image->GetLargestPossibleRegion().GetSize()[1] 
      ) 
    {
      Point<double,ImageDimension> point;   
      for(unsigned int id=0;id<ImageDimension;id++)
      {
        point[id]=y(id);
      }
      v = (m_DataOp->Evaluate(point)-m_DataMin)/(m_DataMax-m_DataMin);
      if(!m_ExtractRidge)
      {
        v = 1-v;
      }
      if(v<0) { v = 0; }
      if(v>1) { v = 1; }
      kernNeg[j] = v*v;
    }
    j++;
  }
    
  int k, kernNegMinK, kernNegMaxK;
  kernNegMinK = kernNegMaxK = 0;
  double kernNegMin, kernNegMax, kernNegAvg;
  kernNegMax = kernNegMin = kernNegAvg = kernNeg[0];

  double filtTot, kernTot[20];
  kernTot[0] = kernNeg[0]/kernPos;
  filtTot = kernTot[0];    
  for(k=1; k<j; k++) 
  {
    kernTot[k] = kernNeg[k]/kernPos;
    filtTot += kernTot[k];
    kernNegAvg += kernNeg[k];
    if(kernNeg[k]>kernNegMax) 
    {
      kernNegMax = kernNeg[k];
      kernNegMaxK = k;       
    }
    if(kernNeg[k]<kernNegMin) 
    {
      kernNegMin = kernNeg[k];
      kernNegMinK = k;
    } 
  }
   
  k = 1;
  filtTot -= kernTot[kernNegMinK];
  kernNegAvg -= kernNegMin;
    
  if(kernNegMaxK != kernNegMinK) 
  {
    filtTot -= kernTot[kernNegMaxK];
    kernNegAvg -= kernNegMax;
    k++;
  }
  filtTot /= (j-k);
  kernNegAvg /= (j-k);
    
  // pnt->SetBranchness((kernNegMax/kernPos)/(kernNegAvg/kernPos));   
  pnt.SetBranchness((kernNegMax/kernPos)/(kernNegAvg/kernPos));    
  return (kernNegMax/kernPos)/(kernNegAvg/kernPos);

}

/**
 * Compute the medialness at a kernel */
template<class TInputImage>
double 
RadiusExtractor2D<TInputImage>
::MedialnessAtKern(std::list<TubePointType>  * tube, double pntR)
{
  if(tube->size() == 0)
  {
    return 0;
  }
   
  if(pntR<m_MedialnessOpt.xMin())
  {
    pntR = m_MedialnessOpt.xMin();
  } 
  if(pntR>m_MedialnessOpt.xMax())
  {
    pntR = m_MedialnessOpt.xMax();
  } 
    
  int i;
  int len = tube->size();
  double w, wTot = 0;
  m_KernMedial = 0;
    
  if(m_Debug) 
  {
    std::cout << "RadiusExtractor2D: medialnessAtCurve: pntR = " << pntR << " : size = " << tube->size() << std::endl;
  }

  std::list<TubePointType>::iterator pnt = tube->begin();
  
  for(i=0, pnt = tube->begin(); pnt != tube->end(); pnt++, i++)
  {
    w = 1.0-fabs(i-len/2.0)/(len/2.0+1);

    double TempMedialness =  MedialnessAtPoint(*pnt, pntR);
    m_KernMedial += w * TempMedialness;
    wTot += w;
    if(m_Debug) 
      {
       std::cout << " i = " << i << " : w = " << w ;
      std::cout << "   M = " << TempMedialness << std::endl;
      }
  }
    
  m_KernMedial /= wTot;

  if(m_Debug)
  {
    std::cout << " totM = " << m_KernMedial << std::endl;
  }
    
  return m_KernMedial;

}

/**
 * Compute the Branchness at a kernel */
template<class TInputImage>
double
RadiusExtractor2D<TInputImage> 
::BranchnessAtKern(TubeType * tube, double pntR)
{
  if(tube->GetPoints()->size() == 0)
  {
    return 0;
  }
   
  if(pntR<m_MedialnessOpt.xMin())
  {
    pntR = m_MedialnessOpt.xMin();
  }
  if(pntR>m_MedialnessOpt.xMax())
  {
    pntR = m_MedialnessOpt.xMax();
  } 
    
  int i;
  int len = tube->GetPoints()->size();
  double w, wTot = 0;
  m_KernBranch = 0;
    
  if(m_Debug) 
  {
    std::cout << "RadiusExtractor2D: medialnessAtCurve: pntR = " << pntR << " : size = " << tube->GetPoints()->size() << std::endl;
  }
  
  std::list<TubePoint*>::iterator pnt = tube->GetPoints()->begin();
    
  for(i=0, pnt = tube->GetPoints()->begin(); pnt != tube->GetPoints()->end(); pnt++, i++) 
  {
    w = 1.0-fabs(i-len/2.0)/(len/2.0+1);
    double TempBranchness  = BranchnessAtPoint((*pnt), pntR);
   m_KernBranch += w * TempBranchness;
    wTot += w;
    if(m_Debug) 
      {
      std::cout << " i = " << i << std::endl;
      std::cout << "   B = " << TempBranchness << std::endl;
      }
  }
    
  m_KernBranch /= wTot;

  if(m_Debug) 
  {
    std::cout << " totB = " << m_KernBranch << std::endl;
  }

  return m_KernBranch;

}


/**
 * Calculate the Optimal scale
 */
template<class TInputImage>
bool
RadiusExtractor2D<TInputImage>
::CalcOptimalScale(TubePointType pnt, bool itkNotUsed( firstGuess ))
{
  PointType x;
  ITKVectorType t;
   
  m_Kern->clear();
  TubePointType tmpPnt;
    
  x[0] = (pnt.GetPosition())[0] + (pnt.GetTangent())[0];
  x[1] = (pnt.GetPosition())[1] + (pnt.GetTangent())[1];
   
  tmpPnt.SetPosition(x);

  t = pnt.GetTangent();

  tmpPnt.SetTangent(t);
  m_Kern->push_front(tmpPnt);

  tmpPnt.SetPosition(pnt.GetPosition());
  tmpPnt.SetTangent(pnt.GetTangent());
  m_Kern->push_front(tmpPnt);
  x[0] = (pnt.GetPosition())[0] - (pnt.GetTangent())[0];
  x[1] = (pnt.GetPosition())[1] - (pnt.GetTangent())[1];
  tmpPnt.SetPosition(x);
  t = pnt.GetTangent();
  tmpPnt.SetTangent(t);
  m_Kern->push_front(tmpPnt);
  
  double xMin = (m_Radius0-10);
  if(xMin < 0.5)
  {
    xMin = 0.5;
  }
  double xMax = (m_Radius0+10);
  if(xMax > m_MedialnessOpt.xMax())
  {
    xMax = m_MedialnessOpt.xMax();
  }
  double mMaxX = m_Radius0;

  double mMax = MedialnessAtKern(m_Kern, mMaxX);
 
  double dir = 0.5;
  double tf = MedialnessAtKern(m_Kern, mMaxX+dir);
  
  
  if(tf<mMax) 
  {
    dir = -dir;
    tf = MedialnessAtKern(m_Kern, mMaxX+dir);
  }

  while(tf>=mMax && mMaxX!=xMin && mMaxX!=xMax) 
  {
    mMax = tf;
    mMaxX = mMaxX+dir;
    tf = MedialnessAtKern(m_Kern, mMaxX+dir);
  }
  if(mMax<m_ThreshWValStart) 
  {
    pnt.SetRadius(m_Radius0);
    return false;
  }
  
  pnt.SetRadius(mMaxX);
  if( m_Debug)
    {
    std::cout<<"RadiusExtractor2D::CalcOptimalScale: the scale chosen: "
             <<mMaxX<<std::endl;
    }

  return true;
}


/**
 * Compute Radii one way */
template<class TInputImage>
bool
RadiusExtractor2D<TInputImage>
::CalcRadiiOneWay(std::vector<TubePointType>::iterator tubeFromPnt, std::vector<TubePointType>::iterator tubeToPnt, bool forward)
{
  int i, j;
  double pntR = m_Radius0;
  double prevPntR = m_Radius0;
  std::vector<TubePointType>::iterator pnt;
  std::vector<TubePointType>::iterator pntStart, pntEnd;
  pntStart = pntEnd = tubeFromPnt;
  bool newPnt;
  int countPrev = 0;
  double wTot, w;
  bool done = false;
  
  pnt = tubeFromPnt;
  while(!done && pnt != tubeToPnt) 
    {
    m_Kern->clear();
    pntStart = pnt;
    newPnt = false;
    for(i=0; !done && pnt != tubeToPnt && i<m_NumRadiusPoints; i++)
      {
      TubePointType cPnt;
      if(pnt == pntEnd)
        {
        newPnt = true;
        }
      
      if(newPnt)
        {
        m_TubePointCount++;
        }   
      wTot = 0;

      for(j=0; pnt != tubeToPnt && j<m_RadiusPointSpacing; j++) 
        {
        w = 1-fabs(j-m_RadiusPointSpacing/2.0)/(m_RadiusPointSpacing/2.0+1);//it looks like a Guassian

        cPnt.SetPosition((cPnt.GetPosition())[0] + w * ((*pnt).GetPosition())[0],
                         (cPnt.GetPosition())[1] + w * ((*pnt).GetPosition())[1]);
        cPnt.SetTangent( (cPnt.GetTangent())[0] + w * ((*pnt).GetTangent())[0],
                          (cPnt.GetTangent())[1] + w * ((*pnt).GetTangent())[1]);

        wTot += w;
         
        if(m_Debug)
          {
          std::cout << "   " << i << ":" << w << ":" << (cPnt.GetPosition())[0]/wTot << ", " << (cPnt.GetPosition())[1]/wTot<< std::endl;
          }
        if(forward)
          {
          pnt++;
          }
        else
          {
          pnt--;
          }
        }
      
      if(j != m_RadiusPointSpacing)
        {
        done = true;
        }
      
      cPnt.SetPosition((cPnt.GetPosition())[0] / wTot,
                       (cPnt.GetPosition())[1] / wTot);
      cPnt.SetTangent((cPnt.GetTangent())[0] / wTot,
                        (cPnt.GetTangent())[1] / wTot);
      
      ITKVectorType t = cPnt.GetTangent();
      t.Normalize();
      cPnt.SetTangent(t);

      if(m_Debug) 
        {
        std::cout << i << ":" << cPnt.GetPosition()[0] << ", " << cPnt.GetPosition()[1]  << std::endl;
        }
      
      m_Kern->push_front(cPnt);
      }
    
    if(i != m_NumRadiusPoints)
      {
      done = true;
      }

    pntEnd = pnt;
    
    if(m_Debug) 
      {   
      std::cout << "pntR0 = " << pntR << std::endl;
      std::cout << "   len = " << m_Kern->size() << std::endl;
      }

    if(m_MedialnessOpt.extreme(&pntR, &w)==false)
      {
        if(m_Debug) 
          {
          std::cout << " m_MedialnessOpt.extreme false " << std::endl;
          }
      pntR = prevPntR;
      }
  
    if(m_Debug) 
      {
      std::cout << "   pntRFinal = " << pntR << " : w = " << w << std::endl;
      std::cout << "   m_ThreshWVal = " << m_ThreshWVal << std::endl;
      }

    if(w<m_ThreshWVal)
      {
      pntR = prevPntR;
      m_MedialnessOpt.extreme(&pntR, &w);
      if(w>m_ThreshWVal) 
        {
        if(w>2*m_ThreshWVal)
          {
          prevPntR = pntR;
          }
        else
          {
          std::cout << "Kurt : what if here! "<<std::endl;
          }
        }
      else 
        {
        pntR = prevPntR;
        w = (double)(2*m_ThreshWVal);
        }
      }
    else
      {        
      if(w>2*m_ThreshWVal)
        {
        prevPntR = pntR;
        }
      }
        if(m_Debug) 
      {
      std::cout << "   pntRFinalFinal = " << pntR << " : w = " << w << std::endl;
      }
    //BranchnessAtKern(m_Kern, pntR);
    double oldW = w;
     m_KernBranch = 0.0;
    i = m_Kern->size() * m_RadiusPointSpacing;
      
    for(j=0, pnt=pntStart; pnt!=pntEnd; j++) 
      {
      w = (double)(1-(1+fabs(j-i/2.0))/(i/2.0+2));
      if((*pnt).GetRadius()!=0)
        {
        (*pnt).SetRadius((1-w)*(*pnt).GetRadius()+w*pntR);
        }
      else
        {
        if(oldW>m_ThreshWVal)
          {
          (*pnt).SetRadius(pntR);
          }
        }
      if((*pnt).GetMedialness()!=0)
        {
        (*pnt).SetMedialness((1-w)*(*pnt).GetMedialness()+w*m_KernMedial);
       
        }
      else
        {
        (*pnt).SetMedialness(m_KernMedial);
        
        }
      if((*pnt).GetBranchness()<m_KernBranch)
        { 
        (*pnt).SetBranchness(m_KernBranch);
        }
      if(forward)
        {
        pnt++;
        }
      else
        {
        pnt--;
        }
      }

    for(j=0, pnt=pntStart;!done && j<m_NumRadiusPoints && pnt!=pntEnd; j++, m_TubePointCount++)
      {     
      if(forward)
        {
        pnt++;
        }
      else
        {
        pnt--;
        }
      }       
    
    if(fabs((float)m_TubePointCount-countPrev)>50)
      {
      countPrev = m_TubePointCount;
      if(m_IdleCallBack!=NULL)
        {
        m_IdleCallBack();
        }
      if(m_StatusCallBack) 
        {
        char st[80];
        sprintf(st, "%d of %d", m_TubePointCount, m_TubeLength);
        m_StatusCallBack("Extract:Widths", st, 1);
        }
      else
        {
        std::cout << m_TubePointCount << " of " << m_TubeLength << std::endl;
        }
      }

    if(forward)
      {
      pnt++;
      }
    else
      {
      pnt--;
      }
    }

  if(m_StatusCallBack) 
    {
    char s[80];
    sprintf(s, "%d of %d", m_TubeLength, m_TubeLength);
    m_StatusCallBack("Extract:Widths", s, 0);
    }
  else
    {
    std::cout << m_TubeLength << " of " << m_TubeLength << std::endl;
    }
 
  return true;
}


/**
 * Calculate Radii */
template<class TInputImage>
bool
RadiusExtractor2D<TInputImage>
::CalcRadii(TubeType * tube)
{
  if(tube == NULL || tube->GetPoints().size() == 0)
  {
    return false;
  } 
    
  m_TubePointCount = 0;
  m_TubeLength = tube->GetPoints().size();

  std::vector<TubePointType>::iterator pnt;  
  for(pnt = tube->GetPoints().begin(); (*pnt).GetID() != 0 && pnt != tube->GetPoints().end(); pnt++);

  TubePointType tp;
  tp = *pnt;

  std::cout << "Found point " << (*pnt).GetID() << std::endl;

  int i;
  for(i=0; i<m_RadiusPointSpacing*3 && pnt!= tube->GetPoints().begin(); i++, pnt--);
  std::cout << "Starting at point " << (*pnt).GetID() << std::endl;

  CalcRadiiOneWay(pnt, tube->GetPoints().end(), true);
  m_TubePointCount -= m_RadiusPointSpacing*6;

  for(i=0; i<m_RadiusPointSpacing*3 && pnt!= tube->GetPoints().end(); i++, pnt++);
  if(pnt == tube->GetPoints().end())
  {
    pnt--;
  }
  std::cout << "Ending (reverse) at point " << (*pnt).GetID() << std::endl;
   
  CalcRadiiOneWay(pnt, tube->GetPoints().begin(), false);
  std::cout << "point id=0 r=" << tp.GetRadius() << std::endl;
  

  pnt = tube->GetPoints().begin();
  pnt++;

  (*tube->GetPoints().begin()).SetRadius((*pnt).GetRadius());
  (*tube->GetPoints().begin()).SetMedialness((*pnt).GetMedialness());
   
  (*tube->GetPoints().begin()).SetBranchness((*pnt).GetBranchness());
  
  return true;
}

/**
 * Idle callback */
template<class TInputImage>
void
RadiusExtractor2D<TInputImage>
::SetIdleCallBack(bool (*idleCallBack)())
{
  m_IdleCallBack = idleCallBack;
}

/**
 * Satus Call back */
template<class TInputImage>
void
RadiusExtractor2D<TInputImage>
::SetStatusCallBack(void (*statusCallBack)(char *, char *, int))
{
  m_StatusCallBack = statusCallBack;
}


}; // end namespace itk

#endif
