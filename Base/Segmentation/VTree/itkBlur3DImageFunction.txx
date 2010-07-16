/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkBlur3DImageFunction.txx,v $
  Language:  C++
  Date:      $Date: 2005/09/06 21:26:16 $
  Version:   $Revision: 1.9 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkBlur3DImageFunction_txx
#define _itkBlur3DImageFunction_txx

#include "itkBlur3DImageFunction.h"

#include <cmath>
#include <algorithm>

#include "itkImage.h"
#include "itkContinuousIndex.h"


namespace itk
{
/**
 * Set the input Image */
template <class TInputImage>
Blur3DImageFunction<TInputImage>
::Blur3DImageFunction( )
{
  m_Debug = false;

  this->m_Image = NULL;
  m_ImageSize.Fill(0);
  m_Spacing.Fill(0);
  m_OriginalSpacing.Fill(0);
  m_UseRelativeSpacing = true;
  m_Scale = 1;
  m_Extent = 3.1;
}

/**
 * Set the input Image */
template <class TInputImage>
void
Blur3DImageFunction<TInputImage>
::SetInputImage( const InputImageType * ptr )
{
  this->Superclass::SetInputImage( ptr );
  m_ImageSize = 
     this->GetInputImage()->GetLargestPossibleRegion().GetSize();
  m_Spacing  = this->GetInputImage()->GetSpacing();
  m_OriginalSpacing  = this->GetInputImage()->GetSpacing();
  if(m_UseRelativeSpacing)
    {
    for(int i=0; i<ImageDimension; i++)
      {
      m_Spacing[i] = m_OriginalSpacing[i] / m_OriginalSpacing[0];
      }
    }

  /* Values by default */
  this->RecomputeKernel();

}


/**
 * Print */
template<class TInputImage>
void
Blur3DImageFunction<TInputImage>
::SetUseRelativeSpacing(bool useRelativeSpacing) 
{
  m_UseRelativeSpacing = useRelativeSpacing;
  if(m_UseRelativeSpacing)
    {
    for(int i=0; i<ImageDimension; i++)
      {
      m_Spacing[i] = m_OriginalSpacing[i] / m_OriginalSpacing[0];
      }
    }
  else
    {
    for(int i=0; i<ImageDimension; i++)
      {
      m_Spacing[i] = m_OriginalSpacing[i];
      }
    }

}

/**
 * Print */
template<class TInputImage>
void
Blur3DImageFunction<TInputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "calculate Blurring value at point:" << std::endl;
}


/**
 * SetScale
 * Pre-compute kernel weights */
template<class TInputImage>
void
Blur3DImageFunction<TInputImage>
::RecomputeKernel( void )
{
  if(m_Debug)
    {
    std::cout << "RecomputeKernel" << std::endl;
    }
  double gfact = -0.5/(m_Scale*m_Scale);

  for(int i=0; i<ImageDimension; i++)
    {
    m_KernelMax[i] = (int)((m_Scale*m_Extent)/m_Spacing[i]);
    if(m_KernelMax[i]<1)
      {
      m_KernelMax[i] = 1;
      }
    m_KernelMin[i] = -m_KernelMax[i];
    }
  if(m_Debug)
    {
    std::cout << "  Scale = " << m_Scale << std::endl;
    std::cout << "  Extent = " << m_Extent << std::endl;
    std::cout << "  KernelMin = " << m_KernelMin << std::endl;
    std::cout << "  KernelMax = " << m_KernelMax << std::endl;
    }
  
  m_KernelWeights.clear();
  m_KernelX.clear();

  IndexType index;
  m_KernelTotal = 0;
  for(index[2] = m_KernelMin[2]; index[2]<=m_KernelMax[2]; index[2]++)
    {
    double distZ = index[2] * m_Spacing[2];
    distZ = distZ * distZ;
    for(index[1] = m_KernelMin[1]; index[1]<=m_KernelMax[1]; index[1]++)
      {
      double distY = index[1] * m_Spacing[1];
      distY = distY * distY + distZ;
      for(index[0] = m_KernelMin[0]; index[0]<=m_KernelMax[0]; index[0]++)
        {
        double dist = index[0] * m_Spacing[0];
        dist = dist * dist + distY;
        double w = exp(gfact*(dist));
        m_KernelWeights.push_back(w);
        m_KernelX.push_back(index);
        m_KernelTotal += w;
        }
      }
    }
}

/**
 * SetScale
 * Pre-compute kernel weights */
template<class TInputImage>
void
Blur3DImageFunction<TInputImage>
::SetScale(double scale)
{
  if(m_Scale != scale)
    {
    m_Scale = scale;
    this->RecomputeKernel();
    }
}

/**
 * SetExtent
 * Pre-compute kernel weights */
template<class TInputImage>
void
Blur3DImageFunction<TInputImage>
::SetExtent(double extent)
{
  if(m_Extent != extent)
    {
    m_Extent = extent;
    this->RecomputeKernel();
    }
}

/**
 * Evaluate the fonction at the specified point */
template <class TInputImage>
double
Blur3DImageFunction<TInputImage>
::Evaluate(const PointType& point) const
{
  if(m_Debug)
    {
    std::cout << "Blur3DImageFunction::Evaluate" << std::endl;
    }

  ContinuousIndexType index;
  if( !this->m_Image )
    {
    index[0] = point[0];
    index[1] = point[1];
    index[2] = point[2];
    }
  else
    {
    this->m_Image->TransformPhysicalPointToContinuousIndex(point,index);
    }

  if(m_Debug)
    {
    std::cout << "  Calling EvaluateAtContinuousIndex " << std::endl;
    }
  return this->EvaluateAtContinuousIndex(index);
}

template <class TInputImage>
double
Blur3DImageFunction<TInputImage>
::EvaluateAtIndex(const IndexType & point) const
{
  if(m_Debug)
    {
    std::cout << "Blur3DImageFunction::EvaluateAtIndex" << std::endl;
    std::cout << "  Point = " << point << std::endl;
    }

  if( !this->m_Image )
    {
    return 0.0;
    }
  
  double res = 0;
  double wTotal;

  IndexType kernelX;

  bool boundary = false;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    if(point[i]+m_KernelMin[i]<0 || point[i]+m_KernelMax[i]>=(int)m_ImageSize[i])
      {
      boundary = true;
      }
    }

  if(!boundary)
    {
    KernelWeightsListType::const_iterator it;
    KernelWeightsListType::const_iterator itEnd;
    KernelXListType::const_iterator  itX;
    it = m_KernelWeights.begin();
    itEnd = m_KernelWeights.end();
    itX = m_KernelX.begin();
    while(it != itEnd)
      {
      kernelX[0] = point[0] + (*itX)[0];
      kernelX[1] = point[1] + (*itX)[1];
      kernelX[2] = point[2] + (*itX)[2];
      res +=  this->m_Image->GetPixel( kernelX ) * (*it);
      ++it;
      ++itX;
      }
    wTotal = m_KernelTotal; 
    }
  else
    {
    if(m_Debug)
      {
      std::cout << "  Boundary point" << std::endl;
      }
    KernelWeightsListType::const_iterator it;
    KernelWeightsListType::const_iterator itEnd;
    KernelXListType::const_iterator  itX;
    it = m_KernelWeights.begin();
    itEnd = m_KernelWeights.end();
    itX = m_KernelX.begin();
    wTotal = 0; 
    double w;
    while(it != itEnd)
      {
      kernelX[0] = point[0] + (*itX)[0];
      kernelX[1] = point[1] + (*itX)[1];
      kernelX[2] = point[2] + (*itX)[2];
      if(kernelX[0]>=0 && kernelX[0]<(int)m_ImageSize[0] &&
         kernelX[1]>=0 && kernelX[1]<(int)m_ImageSize[1] &&
         kernelX[2]>=0 && kernelX[2]<(int)m_ImageSize[2])
        {
        w = *it;
        res +=  this->m_Image->GetPixel( kernelX ) * w;
        wTotal += w;
        }
      ++it;
      ++itX;
      }
    }
  
  if(wTotal < 1e-6)
    { 
    return 0;
    }

  if(m_Debug)
    {
    std::cout << "  result = " << res/wTotal << std::endl;
    }
  return res/wTotal;
}

template <class TInputImage>
double
Blur3DImageFunction<TInputImage>
::EvaluateAtContinuousIndex(const ContinuousIndexType & point) const
{
  if(m_Debug)
    {
    std::cout << "Blur3DImageFunction::EvaluateAtContinuousIndex" << std::endl;
    std::cout << "  Point = " << point << std::endl;
    }

  if( !this->m_Image )
    {
    return 0.0;
    }
  
  double w;
  double res = 0;
  double wTotal = 0;
  double gfact = -0.5/(m_Scale*m_Scale);
  //double kernrad = m_Scale*m_Extent*m_Scale*m_Extent;

  IndexType kernelX;

  bool boundary = false;
  for(unsigned int i=0; i<ImageDimension; i++)
    {
    if(point[i]+m_KernelMin[i]<0 || point[i]+m_KernelMax[i]>=m_ImageSize[i])
      {
      boundary = true;
      }
    }

  if(!boundary)
    {
    for(int z = m_KernelMin[2]; z<=m_KernelMax[2]; z++)
      {
      kernelX[2] = (int)(point[2])+z;
      double distZ = (kernelX[2]-point[2])*m_Spacing[2];
      distZ = distZ * distZ;
      for(int y = m_KernelMin[1]; y<=m_KernelMax[1]; y++)
        {
        kernelX[1] = (int)(point[1])+y;
        double distY = (kernelX[1]-point[1])*m_Spacing[1];
        distY = distY * distY + distZ;
        for(int x = m_KernelMin[0]; x<=m_KernelMax[0]; x++)
          {
          kernelX[0] = (int)(point[0])+x;
          double distX = (kernelX[0]-point[0])*m_Spacing[0];
          double dist = distX * distX + distY;
          w = exp(gfact*dist);
          wTotal += w; 
          res +=  this->m_Image->GetPixel( kernelX ) * w ;     
          }
        }
      }
    }
  else
    {
    if(m_Debug)
      {
      std::cout << "  Boundary point" << std::endl;
      }
    int zero = 0;
    int minZ = vnl_math_min((int)((int)(point[2])+m_KernelMin[2]), zero);
    int minY = vnl_math_min((int)((int)(point[1])+m_KernelMin[1]), zero);
    int minX = vnl_math_min((int)((int)(point[0])+m_KernelMin[0]), zero);
    int maxZ = vnl_math_min((int)((int)(point[2])+m_KernelMax[2]),
                            (int)(m_ImageSize[2]-1));
    int maxY = vnl_math_min((int)((int)(point[1])+m_KernelMax[1]),
                            (int)(m_ImageSize[1]-1));
    int maxX = vnl_math_min((int)((int)(point[0])+m_KernelMax[0]),
                            (int)(m_ImageSize[0]-1));
    for(kernelX[2] = minZ; kernelX[2]<=maxZ; kernelX[2]++)
      {
      double distZ = (kernelX[2]-point[2])*m_Spacing[2];
      distZ = distZ * distZ;
      for(kernelX[1] = minY; kernelX[1]<=maxY; kernelX[1]++)
        {
        double distY = (kernelX[1]-point[1])*m_Spacing[1];
        distY = distY * distY + distZ;
        for(kernelX[0] = minX; kernelX[0]<=maxX; kernelX[0]++)
          {
          double distX = (kernelX[0]-point[0])*m_Spacing[0];
          double dist = distX * distX + distY;
          w = exp(gfact*(dist));
          wTotal += w; 
          res +=  this->m_Image->GetPixel( kernelX ) * w ;     
          }
        }
      }
    }
  
  if(wTotal < 1e-6)
    { 
    return 0;
    }
  if(m_Debug)
    {
    std::cout << "  result = " << res/wTotal << std::endl;
    }
  return res/wTotal;
}

/**
 * Evaluate the fonction at the specified point
 * /warning Need to be rewritten using point dim and point type
 */
template <class TInputImage>
double
Blur3DImageFunction<TInputImage>
::EvaluateLinearInterpolate(const ContinuousIndexType & point) const
{ 
  if(m_Debug)
    {
    std::cout << "Blur3DImageFunction::EvaluateLinearInterpolate" << std::endl;
    std::cout << "  Point = " << point << std::endl;
    }

  int xl1 = (int)point[0];
  int xl2 = (int)point[1];
  int xl3 = (int)point[2];
  double xr1 = point[0] - xl1;
  double xr2 = point[1] - xl2;
  double xr3 = point[2] - xl3;
  int xh1 = xl1+1;
  int xh2 = xl2+1;
  int xh3 = xl3+1;
  if(xh1 >= m_ImageSize[0])   
    {
    xh1 = m_ImageSize[0]-1;
    }
  if(xh2 >= m_ImageSize[1])   
    {
    xh2 = m_ImageSize[1]-1;
    }
  if(xh3 >= m_ImageSize[2]) 
    {
    xh3 = m_ImageSize[2]-1;
    }

  IndexType index1;
  IndexType index2;
  index1[0]=xl1;
  index1[1]=xl2;
  index1[2]=xl3;
  
  index2[0]=xl1;
  index2[1]=xl2;
  index2[2]=xh3;

  double vll3 = (1-xr3)* this->m_Image->GetPixel(index1) + xr3*this->m_Image->GetPixel(index2);

  index1[0]=xl1;
  index1[1]=xh2;
  index1[2]=xl3;
  
  index2[0]=xl1;
  index2[1]=xh2;
  index2[2]=xh3;

  double vlh3 = (1-xr3)* this->m_Image->GetPixel(index1) + xr3*this->m_Image->GetPixel(index2);

  index1[0]=xh1;
  index1[1]=xl2;
  index1[2]=xl3;
  
  index2[0]=xh1;
  index2[1]=xl2;
  index2[2]=xh3;

  double vhl3 = (1-xr3)* this->m_Image->GetPixel(index1) + xr3*this->m_Image->GetPixel(index2);

  index1[0]=xh1;
  index1[1]=xh2;
  index1[2]=xl3;
  
  index2[0]=xh1;
  index2[1]=xh2;
  index2[2]=xh3;

  double vhh3 = (1-xr3)* this->m_Image->GetPixel(index1) + xr3*this->m_Image->GetPixel(index2);

  double vl2 = (1-xr2)*vll3 + xr2*vlh3;
  double vh2 = (1-xr2)*vhl3 + xr2*vhh3;
  double v1 = (1-xr1)*vl2 + xr1*vh2;

  if(m_Debug)
    {
    std::cout << "  result = " << v1 << std::endl;
    }

  return v1;
}


} // namespace itk

#endif
