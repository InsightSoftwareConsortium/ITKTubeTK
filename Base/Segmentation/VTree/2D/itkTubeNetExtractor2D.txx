/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkTubeNetExtractor2D.txx,v $
  Language:  C++
  Date:      $Date: 2005/09/16 23:36:26 $
  Version:   $Revision: 1.8 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTubeNetExtractor2D_txx
#define __itkTubeNetExtractor2D_txx

#include "itkTubeNetExtractor2D.h"

namespace itk
{
 
/**
 * Constructor */
template<class TInputImage>
TubeNetExtractor2D<TInputImage>
::TubeNetExtractor2D()
{
  m_TubeNum = 0;
  m_AEUseMask = false;
  m_AEThresh = NULL;
  m_TubeNet = TubeNet::New();
}

/**
 * Destructor */
template<class TInputImage>
TubeNetExtractor2D<TInputImage>
::~TubeNetExtractor2D()
{
  if(m_AEThresh)
  {
    delete m_AEThresh;
  }
}

/**
 * Set the input image */
template<class TInputImage>
void
TubeNetExtractor2D<TInputImage>
::SetInputImage(ImagePointer inputImage )
{
  m_Image = inputImage;
  m_TubeNum = 0;

  m_TubeNet->GetTubes()->clear();
  m_TubeNet->SetNumDimensions(3);
 
  vnl_vector<unsigned int> size_(3);
  typename ImageType::SizeType imageSize = m_Image->GetLargestPossibleRegion().GetSize();

  for(unsigned int i=0;i<ImageDimension;i++)
  {
    size_(i) = imageSize[i];
  }
  
  m_TubeNet->SetDimSize(size_);
  
  TubeExtractor<TInputImage>::SetInputImage(m_Image);
  
  if(m_AEThresh)
  {
    delete m_AEThresh;
  }

  unsigned int size = m_Image->GetLargestPossibleRegion().GetSize()[ImageDimension-1];
  m_AEThresh = new float[size];
  for(unsigned int i=0; i<size; i++)
  {
    m_AEThresh[i] = 0;
  }
}

/**
 * Set the status callback  */
template<class TInputImage>
void
TubeNetExtractor2D<TInputImage>
::NewTubeCallBack(void (*newTubeCallBack)(Tube* ))
{
  //m_TubeNet->NewTubeCallBack(newTubeCallBack);
  TubeExtractor<TInputImage>::NewTubeCallBack(newTubeCallBack);
}

/**
 * Extract a 2D tube */
template<class TInputImage>
bool 
TubeNetExtractor2D<TInputImage>
::ExtractTube(float x, float y)
{

  if(TubeExtractor<TInputImage>::ExtractTube(x, y, m_TubeNum)) 
  {
    if(!this->m_Tube)
    {
      std::cout << "Extract Tube: TubeExtractor returns no tube !" << std::endl;
      return false;
    }   
    m_TubeNet->GetTubes()->push_back( this->m_Tube);
    m_TubeNum++;
    return true;
  }
  return false;

}


/**
 * Delete a tube */
template<class TInputImage>
bool 
TubeNetExtractor2D<TInputImage>
::DeleteTube(Tube * newTube)
{
  TubeExtractor<TInputImage>::DeleteTube(newTube);
  return m_TubeNet->DeleteTube(newTube);
}


/**
 * Get the tube net */
template<class TInputImage>
TubeNet::Pointer
TubeNetExtractor2D<TInputImage>
::GetTubeNet(void)
{
  return m_TubeNet;
}

/**
 * Set the mask  */ 
template<class TInputImage>   
void
TubeNetExtractor2D<TInputImage>
::SetAutoExtractMask(ImagePointer autoExtractMask)
{
  m_AEMask = autoExtractMask;
}

/**
 * Auto extract tubes using a mask */    
template<class TInputImage>   
double  
TubeNetExtractor2D<TInputImage>
::AutoExtractThresh(void)
{
  if(!m_AEThresh)
  {
    return 0;
  }
  return m_AEThresh[0];
}
  
/**
 * Auto extract tubes using a mask */      
template<class TInputImage>   
void   
TubeNetExtractor2D<TInputImage>
::AutoExtractThresh(double newAEThresh)
{
 if(!m_AEThresh)
 {
   return;
 }
 for(unsigned int i=0; i<m_Image->GetLargestPossibleRegion().GetSize()[ImageDimension-1]; i++)
 {     
   m_AEThresh[i] = newAEThresh;
 }
}
  
/**
 * Auto extract tubes using a mask */  
template<class TInputImage>      
void   
TubeNetExtractor2D<TInputImage>
::AutoExtractAutoThresh(double alpha)
{
  if(!m_AEThresh)
  {
    return;
  }
   
  vnl_vector<int> bin(256,0.0);
  int i, j, k, l;
  long n;
  int bnd, tot;
  unsigned int size[ImageDimension];
  for(unsigned int i=0;i<ImageDimension;i++)
  {
    size[i]=m_Image->GetLargestPossibleRegion().GetSize()[i];
  }
  for(i=0; i<size[2]; i++)
  {
    bin = 0.0;
    n = 0;
    for(j=0; j<size[1]; j++)
    {
      for(k=0; k<size[0]; k++)
      {
        typename ImageType::IndexType index;
        index[0]=k;
        index[1]=j;
        index[2]=i;
        if(!m_AEUseMask || m_AEMask->GetPixel(index)>0)
        {
          l = (m_Image->GetPixel(index)- this->m_RidgeOp->GetDataMin())
              /(this->m_RidgeOp->GetDataMax()-this->m_RidgeOp->GetDataMin())*255;
         
          if(l<0)
          {
            l = 0;
          }
          if(l>255)
          {
            l = 255;
          } 
          bin[l]++;
          n++;
        }
      }
    }
    bnd = (int)(alpha*n);
    if(n > 0)
    {
      tot=0;
      for(l=255; l>=0; l--) 
      {
        tot += bin[l];
        if(tot>bnd)
        break;
      }
      m_AEThresh[i] = (l+0.5)/256.0*(this->m_RidgeOp->GetDataMax()
                      -this->m_RidgeOp->GetDataMin())+this->m_RidgeOp->GetDataMin();
    }
    else
    {
      m_AEThresh[i] = this->m_RidgeOp->GetDataMax();
    }
    std::cout << "thesh[" << i << "] = " << m_AEThresh[i] << std::endl;
  }
}
  
/**
 * AutoExtract tubes */
template<class TInputImage>   
bool
TubeNetExtractor2D<TInputImage>   
::AutoExtract(int zMin, int zMax)
{
  this->m_RidgeOp->SetThreshEVRatioStart(-0.325);
  unsigned int size[ImageDimension];
  for(unsigned int i=0;i<ImageDimension;i++)
  {
    size[i]=m_Image->GetLargestPossibleRegion().GetSize()[i];
  }
  char s[80];
  int i, j, k;
  for(i=zMin; i<zMax; i++)
  {
    for(j=5; j<size[1]-5; j++)
    {
      for(k=5; k<size[0]-5; k++)
      { 
        bool extract=false;
        typename ImageType::IndexType index;
        index[0]=k;
        index[1]=j;
        index[2]=i;
        if(m_AEMask->GetPixel(index)>0) // k,j,i
        {
          index[0]+=5;
          if(m_AEMask->GetPixel(index)>0) //k+5,j,i
          {
            index[0]-=10;
            if(m_AEMask->GetPixel(index)>0) //k-5,j,i
            {
              index[0]+=5;
              index[1]+=5;
              if(m_AEMask->GetPixel(index)>0) //k,j+5,i
              {
                index[1]-=10;
                if(m_AEMask->GetPixel(index)>0) //k,j-5,i 
                {
                  index[0]+=5;
                  if(m_AEMask->GetPixel(index)>0) //k+5,j-5,i
                  {
                    index[0]-=10;
                    if(m_AEMask->GetPixel(index)>0)//k-5,j-5,i
                    {
                      index[1]+=10;
                      if(m_AEMask->GetPixel(index)>0)//k-5,j+5,i
                      {
                        index[0]+=10;
                        if(m_AEMask->GetPixel(index)>0)//k+5,j+5,i
                        {
                          extract = true;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        
        index[0]=k;
        index[1]=j;
        index[2]=i;
        
        if( (!m_AEUseMask || extract) 
            && m_Image->GetPixel(index)>m_AEThresh[i]
            && this->m_RidgeOp->GetDataMask()->GetPixel(index) == 0 
          )
        {
          ExtractTube(k, j, i);
        }
      }  
      
      if(j/10 == j/10.0 && this->m_StatusCallBack) 
      {
        sprintf(s, "%d.%d of %d.%d", i, j, size[2]-1,size[1]-1);
        this->m_StatusCallBack("Auto Extract", s, 0);
        if(this->m_IdleCallBack)
        { 
          this->m_IdleCallBack();
        }
      }
    }
  }   
  if(this->m_StatusCallBack)
  {
    this->m_StatusCallBack("Auto Extract", "Done!", 1);
  } 
  return true;
}



}; // end namespace itk

#endif
