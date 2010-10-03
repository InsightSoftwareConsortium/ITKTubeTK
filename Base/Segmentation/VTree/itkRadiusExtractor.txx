/*=========================================================================

Library:   TubeTK/VTree3D

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
#ifndef __itkRadiusExtractor_txx
#define __itkRadiusExtractor_txx

#include "itkRadiusExtractor.h"
#include <itkMinimumMaximumImageFilter.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

namespace itk
{

/** Define the Medialness Function */
template <class T>
class RadiusExtractorMedialnessFunc : public UserFunc<double, double> 
{
  public:

    typedef itk::VesselTubeSpatialObject<3>  TubeType;
    typedef typename TubeType::TubePointType TubePointType;

    std::list<TubePointType>               * m_Tube;

    RadiusExtractor<T> * m_RadiusExtractor;

    RadiusExtractorMedialnessFunc(RadiusExtractor<T> * newRadiusExtractor,
                                  std::list<TubePointType> * newTube)
      {
      m_Tube = newTube;
      m_RadiusExtractor = newRadiusExtractor;
      };
    double value(double x)
      {
      return m_RadiusExtractor->MedialnessAtKern(m_Tube, x, false);
      };
};


/** Constructor */
template<class TInputImage>
  RadiusExtractor<TInputImage>
::RadiusExtractor()
{
  m_Debug = false; 
  m_Verbose = true; 

  m_DataOp = Blur3DImageFunction<ImageType>::New();
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

  m_KernN0.Fill(0);
  m_KernN1.Fill(0);
  m_KernN0[0] = 1;
  m_KernN1[1] = 1;

  m_IterPntArray = new std::vector<TubePointType>::iterator[5000]; 
  m_KernPntArray = new TubePointType[5000]; 
  m_ArrayLen = 0;

  m_Kern.clear();

  m_MedialnessAtKern = new RadiusExtractorMedialnessFunc<TInputImage>(this, & m_Kern);
  m_MedialnessOpt.use(m_MedialnessAtKern);
  m_MedialnessOpt.xMin(0.4);
  m_MedialnessOpt.xMax(20.0);

  m_MedialnessOpt.tolerance(0.001); //TMI:0.001 // 0.005 - 0.025
  m_MedialnessOpt.xStep(0.25); // TMI: 0.1 // Org: 0.49; Good: 0.25; NEXT: 1.0

  m_MedialnessOpt.searchForMin(false);

  m_IdleCallBack = NULL;
  m_StatusCallBack = NULL;
}

/** Destructor */
template<class TInputImage>
  RadiusExtractor<TInputImage>
::~RadiusExtractor()
{
  if(m_MedialnessAtKern != NULL)
    {
    delete m_MedialnessAtKern;
    } 
  m_MedialnessAtKern = NULL;

  m_ArrayLen = 0;
  delete [] m_IterPntArray;
  delete [] m_KernPntArray;
}


/** Set the scale factor */
template<class TInputImage>
void
  RadiusExtractor<TInputImage>
::SetScale(double scale)
{
  m_Scale = scale;
  m_DataOp->SetScale(scale);
}


/** Set the extent factor */
template<class TInputImage>
void
  RadiusExtractor<TInputImage>
::SetExtent(double extent)
{
  m_Extent = extent;
  m_DataOp->SetExtent(extent);
}


/** Set Radius Min */
template<class TInputImage>
void
  RadiusExtractor<TInputImage>
::SetRadiusMin(double radiusMin)
{
  m_RadiusMin = radiusMin;
  m_MedialnessOpt.xMin(m_RadiusMin);
}

/** Set Radius Max */
template<class TInputImage>
void
  RadiusExtractor<TInputImage>
::SetRadiusMax(double radiusMax)
{
  this->m_RadiusMax = radiusMax;
  m_MedialnessOpt.xMax(this->m_RadiusMx);
}

/** Get the medialness operator */
template<class TInputImage>
OptParabolicFit1D & 
  RadiusExtractor<TInputImage>
::GetMedialnessOpt(void)
{
  return & m_MedialnessOpt;
}

/** Set the input image */
template<class TInputImage>
void
  RadiusExtractor<TInputImage>
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
      std::cout << "RadiusExtractor: SetInputImage: Minimum = " << m_DataMin << std::endl;
      std::cout << "RadiusExtractor: SetInputImage: Maximum = " << m_DataMax << std::endl;
      }
    }
}

/**
 * Compute the medialness and the branchness */
template<class TInputImage>
void
RadiusExtractor<TInputImage>
::ComputeMnessBness(double pntR, double w,
                    double *kernPos, double *kernPosCnt,
                    double *kernNeg, double *kernNegCnt,
                    double *kernBrn, double *kernBrnCnt,
                    double & mness, double & bness,
                    bool doBNess)
{
  int i;

  int kernAvgCnt = 0;   
  for(i=0; i<2*m_KernNumT; i++)
    {    
    if(kernNegCnt[i]>=w && kernPosCnt[i]>=w)
      {
      kernAvgCnt++;
      }
    }

  if(kernAvgCnt<2)
    {
    mness = 0;
    bness = 0;
    return;
    }

  int kernNegMinI, kernNegMaxI;
  double kernNegMin, kernNegMax;
  double kernPosMin, kernPosMax;
  double kernNegAvg, kernPosAvg;

  for(i=0; i<2*m_KernNumT; i++)
    {
    if(kernNegCnt[i]>=w && kernPosCnt[i]>=w)
      {     
      break;
      }
    }

  kernNegMinI = i;
  kernNegMaxI = i;
  kernNegMin = kernNeg[i];
  kernNegMax = kernNegMin;
  kernNegAvg = kernNegMin;
  kernPosMin = kernPos[i];
  kernPosMax = kernPosMin;
  kernPosAvg = kernPosMin;

  for(i++; i<2*m_KernNumT; i++) 
    {
    if(kernNegCnt[i]>=w && kernPosCnt[i]>=w)
      {
      if(kernNeg[i]>kernNegMax)
        {
        kernNegMax = kernNeg[i];
        kernNegMaxI = i;
        }
      if(kernNeg[i]<kernNegMin)
        {
        kernNegMin = kernNeg[i];
        kernNegMinI = i;
        }
      if(kernPos[i]>kernPosMax)
        {
        kernPosMax = kernPos[i];
        }
      if(kernPos[i]<kernPosMin)
        {
        kernPosMin = kernPos[i];
        }
      kernNegAvg += kernNeg[i];
      kernPosAvg += kernPos[i];
      }
    }

  if(fabs(kernNegMin-kernNegAvg/kernAvgCnt) >
    fabs(kernNegMax-kernNegAvg/kernAvgCnt) && kernAvgCnt>2)
    {
    kernPosAvg -= kernPos[kernNegMinI];
    kernNegAvg -= kernNegMin;

    kernPosAvg += (kernPos[kernNegMinI] + kernPosAvg/(kernAvgCnt-1)) / 2;
    kernNegAvg += (kernNegMin + kernNegAvg/(kernAvgCnt-1)) / 2;
    int l,m;
    l = (kernNegMinI+1)%(2*m_KernNumT);
    m = (kernNegMinI+2*m_KernNumT-1)%(2*m_KernNumT);
    if(kernNegCnt[l]>=w && kernPosCnt[l]>=w 
      && kernNeg[l]<kernNeg[m])
      {
      kernPosAvg -= kernPos[l];
      kernNegAvg -= kernNeg[l];
      kernPosAvg += (kernPos[l] + kernPosAvg/(kernAvgCnt-1)) / 2;
      kernNegAvg += (kernNeg[l] + kernNegAvg/(kernAvgCnt-1)) / 2;
      }
    else if(kernNegCnt[m]>=w && kernPosCnt[m]>=w
      && kernNeg[l]>=kernNeg[m])
      {
      kernPosAvg -= kernPos[m];
      kernNegAvg -= kernNeg[m];
      kernPosAvg += (kernPos[m]+ kernPosAvg/(kernAvgCnt-1)) / 2;
      kernNegAvg += (kernNeg[m]+ kernNegAvg/(kernAvgCnt-1)) / 2;
      }
    }
  else if(kernAvgCnt>2)
    {
    kernPosAvg -= kernPos[kernNegMaxI];
    kernNegAvg -= kernNegMax;
    kernPosAvg += (kernPos[kernNegMaxI] + kernPosAvg/(kernAvgCnt-1)) / 2;
    kernNegAvg += (kernNegMax + kernNegAvg/(kernAvgCnt-1)) / 2;
    int l, m;
    l = (kernNegMaxI+1)%(2*m_KernNumT);
    m = (kernNegMaxI+2*m_KernNumT-1)%(2*m_KernNumT);
    if(kernNegCnt[l]>=w && kernPosCnt[l]>=w 
      && kernNeg[l]>kernNeg[m])
      {
      kernPosAvg -= kernPos[l];
      kernNegAvg -= kernNeg[l];
      kernPosAvg += (kernPos[l]+ kernPosAvg/(kernAvgCnt-1))/2;
      kernNegAvg += (kernNeg[l]+ kernNegAvg/(kernAvgCnt-1))/2;
      }
    else if(kernNegCnt[m]>=w && kernPosCnt[m]>=w
      && kernNeg[l]<=kernNeg[m])
      {
      kernPosAvg -= kernPos[m];
      kernNegAvg -= kernNeg[m];
      kernPosAvg += (kernPos[m]+ kernPosAvg/(kernAvgCnt-1)) / 2;
      kernNegAvg += (kernNeg[m]+ kernNegAvg/(kernAvgCnt-1)) / 2;
      }
    }
  if(kernAvgCnt!=0)
    {
    kernPosAvg /= kernAvgCnt;
    kernNegAvg /= kernAvgCnt;
    }
  else
    {
    kernPosAvg =0;
    kernNegAvg= 0;
    }  
  mness = (kernPosAvg-kernNegAvg) / (pntR/4);

  if(doBNess)
    {
    double kernBrnMax, kernBrnAvg;
    double kernBrnAvgCnt = 1;
    for(i=0; i<2*m_KernNumT; i++)
      if(kernPosCnt[i]>=w && kernNegCnt[i]>=w && kernBrnCnt[i]>=w)
        break;

    if(i<2*m_KernNumT)
      {
      kernBrnMax = kernBrn[i];
      kernBrnAvg = kernBrnMax;

      for(i++; i<2*m_KernNumT; i++)
        {
        if(kernPosCnt[i]>=w && kernNegCnt[i]>=w && kernBrnCnt[i]>=w)
          {
          kernBrnAvg += kernBrn[i];
          kernBrnAvgCnt++;
          if(kernBrn[i]>kernBrnMax)
            {
            kernBrnMax = kernBrn[i];
            }
          }
        }

      kernBrnAvg -= kernBrnMax;
      if((kernBrnAvgCnt-1)!=0)
        {
        kernBrnAvg /= (kernBrnAvgCnt-1);
        }
      else
        {
        kernBrnAvg=0;
        }

      if((kernPosAvg - kernBrnAvg)!=0)
        {  
        bness = (3 * (kernBrnMax - kernBrnAvg) / (kernPosAvg - kernBrnAvg)
          + (kernNegMax - kernNegAvg) / (kernPosAvg - kernBrnAvg)) / 4;
        }
      else
        {
        bness=0;
        }

      if(bness>2)
        {
        bness = 2;
        }
      if(bness<0)
        {
        bness = 0;
        }
      }
    else
      {
      bness = 0;
      }
    }
  else
    {
    bness = 0;
    }
}


/**
 * Compute the medialness at a point 
 */
template<class TInputImage>
double
  RadiusExtractor<TInputImage>
::MedialnessAtPoint(TubePointType pnt, double pntR, bool doBNess, bool newKern, double w)
{
  if(pntR < m_MedialnessOpt.xMin())
    {
    pntR = m_MedialnessOpt.xMin();
    }

  VectorType y(ImageDimension);
  VectorType n0(ImageDimension);
  VectorType n1(ImageDimension);


  n0.SetVnlVector( GetOrthogonalVector( pnt.GetTangent().GetVnlVector()) );
  n0.Normalize();

  n1.SetVnlVector( itk_cross_3d(pnt.GetTangent().GetVnlVector(), n0.GetVnlVector()) );
  n1.Normalize();

  if(newKern)
    {
    m_KernN0[0] = 1;
    m_KernN0[1] = 0;
    m_KernN0[2] = 0;
    m_KernN1[0] = 0;
    m_KernN1[1] = 1;
    m_KernN1[2] = 0;
    for(int i=0; i<2*this->m_KernNumT; i++)
      {
      m_KernPos[i] = 0;
      m_KernNeg[i] = 0;
      m_KernBrn[i] = 0;
      m_KernPosCnt[i] = 0;
      m_KernNegCnt[i] = 0;
      m_KernBrnCnt[i] = 0;
      }
    }

  double kn0ProjN0 = dot_product(m_KernN0.GetVnlVector(), n0.GetVnlVector()); 
  double kn0ProjN1 = dot_product(m_KernN0.GetVnlVector(), n1.GetVnlVector());
  double kn1ProjN0 = dot_product(m_KernN1.GetVnlVector(), n0.GetVnlVector());
  double kn1ProjN1 = dot_product(m_KernN1.GetVnlVector(), n1.GetVnlVector());

  m_KernN0[0] = kn0ProjN0*n0[0] + kn0ProjN1*n1[0];
  m_KernN0[1] = kn0ProjN0*n0[1] + kn0ProjN1*n1[1];
  m_KernN0[2] = kn0ProjN0*n0[2] + kn0ProjN1*n1[2];

  m_KernN0.Normalize();

  n0 = m_KernN0;

  m_KernN1[0] = kn1ProjN0*n0[0] + kn1ProjN1*n1[0];
  m_KernN1[1] = kn1ProjN0*n0[1] + kn1ProjN1*n1[1];
  m_KernN1[2] = kn1ProjN0*n0[2] + kn1ProjN1*n1[2];

  m_KernN1.Normalize();

  n1.SetVnlVector( itk_cross_3d(pnt.GetTangent().Get_vnl_vector(), n0.GetVnlVector()) );

  double tf = dot_product(m_KernN1.GetVnlVector(), n1.GetVnlVector()); 
  if(tf<0)
    {
    n1[0] = -n1[0];
    n1[1] = -n1[1];
    n1[2] = -n1[2];
    m_KernN1[0] = n1[0];
    m_KernN1[1] = n1[1];
    m_KernN1[2] = n1[2];
    }
  else
    {
    m_KernN1 = n1;
    }

  double kernPos[40];
  double kernNeg[40];
  double kernBrn[40];
  double kernPosCnt[40];
  double kernNegCnt[40];
  double kernBrnCnt[40];
  for(int i=0; i<2*m_KernNumT; i++)
    {
    kernPos[i] = 0;    
    kernNeg[i] = 0;
    kernBrn[i] = 0;
    kernPosCnt[i] = 0;
    kernNegCnt[i] = 0;
    kernBrnCnt[i] = 0;
    }

  int numC, posCMax, negCMax;
  numC = 2;
  if(doBNess)
    {
    numC++;
    }
  posCMax = 1;
  negCMax = 2;

  double v, r;
  for(int c=0; c<numC; c++) 
    {
    if(c < posCMax) // inner (positive) circle
      {
      double e = 1.1;
      if((pntR/4)*e<0.71)  // sqrt(2)/2 = 0.7071 approx = 0.71
        {  
        e = 0.71/(pntR/4); //0.71/(pntR/4);
        }
      if((pntR/4)*e>3.1)
        {      
        e = 3.1 / (pntR/4);
        }   
      m_DataOp->SetScale(pntR/4); 
      m_DataOp->SetExtent(e); 
      r = 3*pntR/4;
      }
    else if(c < negCMax) // outer (negative) circle
      {     
      r = 5*pntR/4;
      }
    else // branchness - wider and futher out
      {
      double e = 1;
      if((pntR/3)*e<1.1)  // sqrt(2)/2 = 0.7071 approx = 0.71
        {
        e = 1.1/(pntR/3);
        }
      if((pntR/3)*e>3.1)
        {
        e = 3.1 / (pntR/3);
        }
      m_DataOp->SetScale(pntR/3); // mess with this -  and r
      m_DataOp->SetExtent(e);
      r = pntR * 3;
      }
    for(int i=0; i<m_KernNumT; i++) 
      {
      VectorType tp;
      tp[0] = pnt.GetPosition()[0];
      tp[1] = pnt.GetPosition()[1];
      tp[2] = pnt.GetPosition()[2];
      y.SetVnlVector( ComputeLineStep(tp.GetVnlVector(), r*m_KernSinT[i], n0.GetVnlVector()) );
      y.SetVnlVector( ComputeLineStep(y.GetVnlVector(), r*m_KernCosT[i], n1.GetVnlVector()) );

      if(y[0]>=0 && y[0]<m_Image->GetLargestPossibleRegion().GetSize()[0] &&
        y[1]>=0 && y[1]<m_Image->GetLargestPossibleRegion().GetSize()[1] &&
        y[2]>=0 && y[2]<m_Image->GetLargestPossibleRegion().GetSize()[2]) 
        {
        ContinuousIndex<double, ImageDimension> point;

        for(unsigned int id=0;id<ImageDimension;id++)
          {
          point[id] = y[id];
          }

        v = (m_DataOp->EvaluateAtContinuousIndex(point)-m_DataMin)/(m_DataMax-m_DataMin); 

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
        if(c < posCMax)
          {
          kernPos[i] += w*v;
          kernPosCnt[i] += w;
          m_KernPos[i] += w*v;
          m_KernPosCnt[i] += w;
          }

        else if(c < negCMax)
          {
          kernNeg[i] += w*v;
          kernNegCnt[i] += w;
          m_KernNeg[i] += w*v;
          m_KernNegCnt[i] += w;
          } 
        else
          {
          kernBrn[i] += w*v;
          kernBrnCnt[i] += w;
          m_KernBrn[i] += w*v;
          m_KernBrnCnt[i] += w;
          }
        }

      tp[0] = pnt.GetPosition()[0];
      tp[1] = pnt.GetPosition()[1];
      tp[2] = pnt.GetPosition()[2];
      y.SetVnlVector( ComputeLineStep(tp.GetVnlVector(), -r*m_KernSinT[i], n0.GetVnlVector()) );
      y.SetVnlVector( ComputeLineStep(y.GetVnlVector(), -r*m_KernCosT[i], n1.GetVnlVector()) );
      if(y[0]>=0 && y[0]<m_Image->GetLargestPossibleRegion().GetSize()[0] &&
         y[1]>=0 && y[1]<m_Image->GetLargestPossibleRegion().GetSize()[1] &&
         y[2]>=0 && y[2]<m_Image->GetLargestPossibleRegion().GetSize()[2]) 
        {
        ContinuousIndex<double, ImageDimension> point;

        for(unsigned int id=0;id<ImageDimension;id++)
          {
          point[id] = y[id];
          }

        v = (m_DataOp->EvaluateAtContinuousIndex(point)-m_DataMin)/(m_DataMax-m_DataMin); 

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

        if(c < posCMax)
          {
          kernPos[m_KernNumT+i] += w*v;
          kernPosCnt[m_KernNumT+i] += w;
          m_KernPos[m_KernNumT+i] += w*v;
          m_KernPosCnt[m_KernNumT+i] += w;
          }
        else if(c < negCMax)
          {
          kernNeg[m_KernNumT+i] += w*v;
          kernNegCnt[m_KernNumT+i] += w;
          m_KernNeg[m_KernNumT+i] += w*v;
          m_KernNegCnt[m_KernNumT+i] += w;
          }
        else
          {
          kernBrn[2*m_KernNumT+i] += w*v;
          kernBrnCnt[2*m_KernNumT+i] += w;
          m_KernBrn[2*m_KernNumT+i] += w*v;
          m_KernBrnCnt[2*m_KernNumT+i] += w;
          }
        }
      }
    }

  int kernAvgCnt = 0;
  for(int i=0; i<2*m_KernNumT; i++)
    {
    if(kernNegCnt[i]>0)
      {
      kernNeg[i] /= kernNegCnt[i];
      }
    if(kernPosCnt[i]>0)
      {
      kernPos[i] /= kernPosCnt[i];
      }
    if(kernNegCnt[i]>=w && kernPosCnt[i]>=w)
      {
      kernAvgCnt++;
      }
    }

  if(kernAvgCnt<2)
    {
    std::cout << "Warning: medialness kernel does not intersect image" << std::endl;
    pnt.SetMedialness(0);
    pnt.SetBranchness(0);
    return 0;
    }

  if(doBNess)
    {
    for(int i=0; i<2*m_KernNumT; i++)
      {     
      if(kernBrnCnt[i]>=w)
        {
        kernBrn[i] /= kernBrnCnt[i];
        }
      }
    }

  double mness, bness;

  ComputeMnessBness(pntR, w, kernPos, kernPosCnt,
                             kernNeg, kernNegCnt,
                             kernBrn, kernBrnCnt,
                             mness, bness, doBNess);
  pnt.SetMedialness(mness);
  pnt.SetBranchness(bness);

  return mness;
}

/** Compute the medialness at a kernel */

template<class TInputImage>
double
  RadiusExtractor<TInputImage>
::MedialnessAtKern(std::list<TubePointType> * tube, double pntR, bool doBNess)

{ 
  //int len = tube->GetPoints().size();
  int len = tube->size();
  int mid = (len-1)/2;
  static double w[20];
  static bool first = true;
  static int vLen = 0;

  if(first || len != vLen)
    {
    vLen = len;
    first = false;
    int i;
    double wTot = 0;
    if(len>1)
      {
      for(i=0; i<len; i++)
        {
        w[i] = 1.0-fabs((double)(i-mid))/(2*mid); // was 1.25 why not ?!
        wTot += w[i];
        }
      for(i=0; i<len; i++)
        {
        w[i] /= wTot;
        }
      }
    else
      {
      w[0] = 1;
      wTot = 1;
      }
    }


  // if(tube->GetPoints().size() == 0)
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
  m_KernMedial = 0;
  m_KernBranch = 0;
  if(m_Debug) 
    {   
    std::cout << "RadiusExtractor3D: MedialnessAtKern: pntR = " 
      << pntR << " : size = " << tube->size() 
      << std::endl;
    }

  std::list<TubePointType>::iterator pnt;
  for(i=0, pnt = tube->begin();
    pnt != tube->end();
    pnt++, i++) 
    {
    MedialnessAtPoint(*pnt, pntR, doBNess, (i==0)?true:false, w[i]);
    }

  for(i=0; i<2*m_KernNumT; i++)
    {
    if(m_KernNegCnt[i]>0)
      {
      m_KernNeg[i] /= m_KernNegCnt[i];
      }
    if(m_KernPosCnt[i]>0)
      {
      m_KernPos[i] /= m_KernPosCnt[i];
      }
    if(m_KernBrnCnt[i]>0)
      {
      m_KernBrn[i] /= m_KernBrnCnt[i];
      }
    }
  ComputeMnessBness(pntR, 0.5, m_KernPos, m_KernPosCnt,
                    m_KernNeg, m_KernNegCnt,
                    m_KernBrn, m_KernBrnCnt, 
                    m_KernMedial, m_KernBranch, doBNess);
                
  return m_KernMedial;
}


/** Compute the Optimal scale */
template<class TInputImage>
bool
  RadiusExtractor<TInputImage>
::CalcOptimalScale(TubePointType pnt, bool firstGuess)
{
  TubePointType tmpPnt;
  PointType x; 

  m_Kern.clear();

  x[0] = pnt.GetPosition()[0] + pnt.GetTangent()[0];
  x[1] = pnt.GetPosition()[1] + pnt.GetTangent()[1];
  x[2] = pnt.GetPosition()[2] + pnt.GetTangent()[2];
  tmpPnt.SetPosition(x);
  tmpPnt.SetTangent(pnt.GetTangent());
  m_Kern.push_front(tmpPnt);

  PointType x1 = pnt.GetPosition();
  tmpPnt.SetPosition(x1);
  tmpPnt.SetTangent(pnt.GetTangent());
  m_Kern.push_front(tmpPnt);

  x1[0] = pnt.GetPosition()[0] + pnt.GetTangent()[0];
  x1[1] = pnt.GetPosition()[1] + pnt.GetTangent()[1];
  x1[2] = pnt.GetPosition()[2] + pnt.GetTangent()[2];
  tmpPnt.SetPosition(x1);
  tmpPnt.SetTangent(pnt.GetTangent());
  m_Kern.push_front(tmpPnt);

  double pntR = m_Radius0;
  double w = 0;

  double tempTol = m_MedialnessOpt.tolerance();
  double tempXStep = m_MedialnessOpt.xStep();
  m_MedialnessOpt.tolerance(0.0001); // 0.005 - 0.025 // TMI 0.05

  if(firstGuess)
    {
    m_MedialnessOpt.xStep(1.0); // 0.05 - 0.25 // TMI 0.25
    }
  else
    {
    m_MedialnessOpt.xStep(0.5); // 0.05 - 0.25 // TMI 0.25
    }

  double tempXMin = m_MedialnessOpt.xMin();
  if(firstGuess)
    {
    m_MedialnessOpt.xMin(0.5);
    if(m_Radius0<1.0)
      {
      m_Radius0 = 1.0;
      }
    m_MedialnessOpt.xMax(m_Radius0*10);
    }

  double tempXMax = m_MedialnessOpt.xMax();

  m_MedialnessOpt.extreme(&pntR, &w);

  m_MedialnessOpt.tolerance(tempTol);
  m_MedialnessOpt.xStep(tempXStep);
  m_MedialnessOpt.xMin(tempXMin);
  m_MedialnessOpt.xMax(tempXMax);

  pnt.SetRadius(pntR);

  if(w>m_ThreshWVal)
    {
    return true;
    }
  else
    {
    std::cout << "RadiusExtractor3D: calcOptimalScale: kernel fit insufficient" << std::endl;
    return false;
    }

  return true;
}


/** Compute kernel array */
template<class TInputImage>
void
  RadiusExtractor<TInputImage>
::CalcKernArray(TubeType * tube)
{
  std::vector<TubePointType>::iterator tubeFromPnt = tube->GetPoints().begin();
  std::vector<TubePointType>::iterator tubeToPnt = tube->GetPoints().end();
  tubeToPnt--;

  std::vector<TubePointType>::iterator iterPnt;

  iterPnt = tubeFromPnt;

  int i, j;
  double wTot;
  i = 0;
  while(iterPnt != tubeToPnt)
    {
    wTot = 0;
    TubePointType kernPnt;
    for(j=0; iterPnt != tubeToPnt && j<m_RadiusPointSpacing; j++) 
      {
      PointType tmpPoint;

      for(unsigned int id=0;id<ImageDimension;id++)
        {
        tmpPoint[id] = kernPnt.GetPosition()[id] + (*iterPnt).GetPosition()[id];
        }

      kernPnt.SetPosition(tmpPoint);
      kernPnt.SetTangent(kernPnt.GetTangent()+(*iterPnt).GetTangent());

      wTot++;
      iterPnt++;
      }

    PointType tmpPoint;
    for(unsigned int id=0;id<ImageDimension;id++)
      {
      tmpPoint[id] = kernPnt.GetPosition()[id] / wTot;
      }

    kernPnt.SetPosition(tmpPoint);
    kernPnt.SetTangent(kernPnt.GetTangent() / wTot);

    VectorType tempVect = kernPnt.GetTangent();
    tempVect.Normalize();
    kernPnt.SetTangent(tempVect);

    kernPnt.SetRadius(0);

    if(iterPnt == tubeToPnt)
      {
      if(i>0)
        {
        //delete kernPnt;
        }
      else
        {
        m_KernPntArray[i] = kernPnt;
        m_IterPntArray[i] = iterPnt;
        i++;
        }
      break;
      }

    m_KernPntArray[i] = kernPnt;
    m_IterPntArray[i] = iterPnt;
    i++;
    }
  m_ArrayLen = i;
}

/** Compute Kernel radii one way */
template<class TInputImage>
void
  RadiusExtractor<TInputImage>
::CalcKernRadiiOneWay(int iStart, int iEnd, bool forward)
{
  int i, j;
  double pntR = m_Radius0;
  double prevPntR = m_Radius0;
  double w;
  char s[80];
  static int count = 0;
  int kernMid = (m_NumRadiusPoints-1)/2;

  if(forward)
    {
    count = 0;
    }

  for(i=iStart; (forward && i<=iEnd) || (!forward && i>=iEnd); (forward)?i++:i--)
    {
    m_Kern.clear();

    for(j=i-kernMid; j<=i+kernMid; j++)
      {
      if((forward && j>=iStart && j<=iEnd) ||
        (!forward && j<=iStart && j>=iEnd))
        {
        m_Kern.push_front(m_KernPntArray[j]); // should be push_front
        }
      else
        {
        m_Kern.push_front(m_KernPntArray[i]); // should be push_front
        }
      }

    m_MedialnessOpt.extreme(&pntR, &w);

    if(w<m_ThreshWVal) 
      {
      if(m_Debug)
        {
        std::cout << "Bad wVal("<<pntR<<") = " << w << std::endl;
        }
      pntR = prevPntR;
      m_MedialnessOpt.extreme(&pntR, &w);
      if(w>m_ThreshWVal) 
        {
        if(m_Debug)
          {
          std::cout << "   *** new wVal("<<pntR<<") = " << w << std::endl;
          }
        if(w>2*m_ThreshWVal)
          { 
          prevPntR = pntR;
          }
        }
      else 
        {
        pntR = prevPntR;
        w = (double)(2*m_ThreshWVal);
        if(m_Debug)
          {
          std::cout << "   using old wVal("<<pntR<<") = " << w << std::endl;
          }
        }
      }
    else
      {
      if(w>2*m_ThreshWVal)
        {
        prevPntR = pntR;
        }
      }

    for(j=i-kernMid; j<=i+kernMid; j++)
      {   
      if((forward && j>=iStart && j<=iEnd) ||
        (!forward && j<=iStart && j>=iEnd))
        {
        if(m_KernPntArray[j].GetRadius()>0)
          {
          w = 1-fabs((double)(j-i))/(double)(kernMid+1);
          w = m_KernPntArray[j].GetRadius() + w * (pntR - m_KernPntArray[j].GetRadius());
          m_KernPntArray[j].SetRadius(w);
          }
        else
          {
          m_KernPntArray[j].SetRadius(pntR);
          }
        }
      }

    count++;
    if(count/5 == count/5.0 && m_StatusCallBack)
      {
      sprintf(s, "Radius at %d = %f", count*m_RadiusPointSpacing, pntR);
      m_StatusCallBack(NULL, s, 1);
      }
    }
}

/** Compute Kernel Measures */
template<class TInputImage>
void
  RadiusExtractor<TInputImage>
::CalcKernMeasures(void)
{
  int i, j;

  for(j=0; j<20; j++)
    {
    for(i=0; i<m_ArrayLen-1; i++)
      {
      m_KernPntArray[i].SetRadius((1.0 * m_KernPntArray[i].GetRadius()
                                   + 1.0 * m_KernPntArray[i+1].GetRadius())/2.0);
      }
    for(i=m_ArrayLen-1; i>0; i--)
      {
      m_KernPntArray[i].SetRadius((1.0 * m_KernPntArray[i].GetRadius()
                                   + 1.0 * m_KernPntArray[i-1].GetRadius())/2.0);
      }
    }

  int kernMid = (m_NumRadiusPoints-1)/2;

  for(i=0; i<m_ArrayLen; i++)
    {
    m_Kern.clear();
    for(j=i-kernMid; j<=i+kernMid; j++)
      {  
      if(j>=0 && j<m_ArrayLen)
        {
        m_Kern.push_front(m_KernPntArray[j]);
        }
      else
        {
        m_Kern.push_front(m_KernPntArray[i]);
        }
      }

    MedialnessAtKern(&m_Kern, m_KernPntArray[i].GetRadius(), true);
    m_KernPntArray[i].SetMedialness(m_KernMedial);
    m_KernPntArray[i].SetBranchness(m_KernBranch);
    }

  m_Kern.clear();
}

/**
 * Apply kernel measures */
template<class TInputImage>
void
  RadiusExtractor<TInputImage>
::ApplyKernMeasures(TubeType * tube)
{

  int i, j;

  double r0 = m_KernPntArray[0].GetRadius();
  double m0 = m_KernPntArray[0].GetMedialness();
  double b0 = m_KernPntArray[0].GetBranchness();

  std::vector<TubePointType> & pnts = tube->GetPoints();
  std::vector<TubePointType>::iterator pnt = pnts.begin();

  while(pnt != pnts.end() && pnt!= m_IterPntArray[0])
    {
    (*pnt).SetRadius(r0);
    (*pnt).SetMedialness(m0);
    std::cout<<"Branch= "<<b0<<std::endl;
    (*pnt).SetBranchness(b0);
    pnt++;
    }
  if(pnt == tube->GetPoints().end())
    {
    return;
    }
  double w, r, m, b;

  j = 0;
  double r1 = r0;
  double m1 = m0;
  double b1 = b0;
  while(j<m_ArrayLen-1)
    {
    r0 = r1;
    m0 = m1;
    b0 = b1;
    r1 = m_KernPntArray[j+1].GetRadius();
    m1 = m_KernPntArray[j+1].GetMedialness();
    b1 = m_KernPntArray[j+1].GetBranchness();
    i = 0;
    while(pnt != m_IterPntArray[j+1])
      {
      w = i/(double)m_RadiusPointSpacing;
      r = r0 + w * (r1 - r0);
      m = m0 + w * (m1 - m0);
      b = b0 + w * (b1 - b0);
      (*pnt).SetRadius(r);
      (*pnt).SetMedialness(m);
      (*pnt).SetBranchness(b);
      pnt++;
      i++;
      }
    j++;
    }

  r0 = m_KernPntArray[j].GetRadius();
  m0 = m_KernPntArray[j].GetMedialness();
  b0 = m_KernPntArray[j].GetBranchness();
  while(pnt != tube->GetPoints().end())
    {
    (*pnt).SetRadius(r0);
    (*pnt).SetMedialness(m0);
    (*pnt).SetBranchness(b0);
    pnt++;
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
  for(i=0; i<m_ArrayLen; i++)
    {
    //    delete m_KernPntArray[i];
    //    m_KernPntArray[i] = NULL;
    //    m_IterPntArray[i] = ;
    }
  return;
}

/** Compute Radii */
template<class TInputImage>
bool
  RadiusExtractor<TInputImage>
::CalcRadii(TubeType * tube)
{
  if(tube == NULL || tube->GetPoints().size() == 0)
    {
    return false;
    } 
  m_TubeLength = tube->GetPoints().size();

  CalcKernArray(tube);

  std::vector<TubePointType>::iterator pnt;
  pnt = tube->GetPoints().begin();
  while( (*pnt).GetID() != 0 && pnt != tube->GetPoints().end() )
    {
    pnt++;
    }

  if(m_Debug)
    {
    std::cout << "Found point " << (*pnt).GetID() << std::endl;
    }

  int i;
  double tf;
  int minDistI = 0;
  double minDist = ComputeEuclideanDistance(m_KernPntArray[0].GetPosition(), (*pnt).GetPosition());


  for(i=1; i<m_ArrayLen; i++)
    {
    tf = ComputeEuclideanDistance(m_KernPntArray[i].GetPosition(),(*pnt).GetPosition());
    if(tf<minDist)
      {
      minDist = tf;
      minDistI = i;
      }
    }

  if(m_Debug)
    {
    std::cout << "Found point i = " << minDistI << std::endl;
    }
  i = minDistI-1;
  if(i<0)
    {
    i = 0;
    }

  CalcKernRadiiOneWay(i, m_ArrayLen-1, true);

  i = minDistI+1;
  if(i>m_ArrayLen-1)
    {
    i = m_ArrayLen-1;
    }

  CalcKernRadiiOneWay(i, 0, false);
  //CalcKernRadiiOneWay(0, i, true);
  CalcKernMeasures();
  ApplyKernMeasures(tube);

  return true;
}

/**
 * Idle callback */
template<class TInputImage>
void
  RadiusExtractor<TInputImage>
::SetIdleCallBack(bool (*idleCallBack)())
{
  m_IdleCallBack = idleCallBack;
}

/**
 * Satus Call back */
template<class TInputImage>
void
  RadiusExtractor<TInputImage>
::SetStatusCallBack(void (*statusCallBack)(char *, char *, int))
{
  m_StatusCallBack = statusCallBack;
}


}; // end namespace itk

#endif

