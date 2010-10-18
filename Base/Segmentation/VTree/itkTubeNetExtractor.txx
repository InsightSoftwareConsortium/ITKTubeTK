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
#ifndef __itkTubeNetExtractor_txx
#define __itkTubeNetExtractor_txx

#include "itkTubeNetExtractor.h"

namespace itk
{

/**
* Constructor */
template<class TInputImage, class TInputMask>
TubeNetExtractor<TInputImage, TInputMask>
::TubeNetExtractor()
{
  m_TubeNum = 0;
  m_AEUseMask = false;
  m_AEThresh = 0;
  m_TubeNet = TubeType::New();
}

/**
 * Destructor */
template<class TInputImage, class TInputMask>
TubeNetExtractor<TInputImage, TInputMask>
::~TubeNetExtractor()
{
}

/**
 * Set the input image */
template<class TInputImage, class TInputMask>
void
TubeNetExtractor<TInputImage, TInputMask>
::SetInputImage(ImagePointer inputImage )
{
  m_Image = inputImage;
  m_TubeNum = 0;

  //m_TubeNet->GetTubes()->clear();
  //m_TubeNet->SetNumDimensions(3);
  /* 
     vnl_vector<unsigned int> size_(3);
     ImageType::SizeType imageSize = m_Image->GetLargestPossibleRegion().GetSize();

     for(unsigned int i=0;i<ImageDimension;i++)
     {
     size_(i) = imageSize[i];
     }

     m_TubeNet->SetDimSize(size_);
   */  
  TubeExtractor<TInputImage>::SetInputImage(m_Image);


  m_AEThresh = 0;
}

/**
 * Set the status callback  */
template<class TInputImage, class TInputMask>
void
TubeNetExtractor<TInputImage, TInputMask>
::NewTubeCallBack(void (*newTubeCallBack)(TubeType* ))
{
  //m_TubeNet.newTubeCallBack(newNewTubeCallBack);
  TubeExtractor<TInputImage>::NewTubeCallBack(newTubeCallBack);
}

/**
 * Extract a 3D tube */
template<class TInputImage, class TInputMask>
bool 
TubeNetExtractor<TInputImage, TInputMask>
::ExtractTube(float x, float y, float z)
{
  /*if(TubeExtractor<TInputImage>::ExtractTube(x, y, z, m_TubeNum)) 
    {
    if(!m_Tube)
      {
      std::cout << "Extract Tube: TubeExtractor returns no tube !" << std::endl;
      return false;
      }   
    m_TubeNet->AddSpatialObject(m_Tube);
    m_TubeNum++;
    return true;
    }*/
  return false;
}

/**
 * Delete a tube */
template<class TInputImage, class TInputMask>
bool 
TubeNetExtractor<TInputImage, TInputMask>
::DeleteTube(TubeType * newTube)
{
  //TubeExtractor<TInputImage>::DeleteTube(newTube);
  //return m_TubeNet->DeleteTube(newTube);
  return false;
}


/**
 * Get the tube net */
template<class TInputImage, class TInputMask>
typename TubeNetExtractor<TInputImage, TInputMask>::TubePointer
TubeNetExtractor<TInputImage, TInputMask>
::GetTubeNet(void)
{
  return m_TubeNet;
}

/**
 * Set the mask  */ 
template<class TInputImage, class TInputMask>   
void
TubeNetExtractor<TInputImage, TInputMask>
::SetAutoExtractMask(MaskPointer autoExtractMask)
{
  m_AEMask = autoExtractMask;
}

/**
 * Auto extract tubes using a mask */    
template<class TInputImage, class TInputMask>   
double  
TubeNetExtractor<TInputImage, TInputMask>
::AutoExtractThresh(void)
{
  return m_AEThresh;
}

/**
 * Auto extract tubes using a mask */      
template<class TInputImage, class TInputMask>   
void   
TubeNetExtractor<TInputImage, TInputMask>
::AutoExtractThresh(double newAEThresh)
{
  m_AEThresh = newAEThresh;
}

/**
 * Auto extract tubes using a mask */  
template<class TInputImage, class TInputMask>      
void   
TubeNetExtractor<TInputImage, TInputMask>
::AutoExtractAutoThresh(double alpha)
{
/*  vnl_vector<int> bin(256,0.0);
  int i, j, k, l;
  long n;
  int bnd, tot;
  unsigned int size[ImageDimension];
  for(unsigned int i=0;i<ImageDimension;i++)
    {
    size[i]=m_Image->GetLargestPossibleRegion().GetSize()[i];
    }
  bin = 0.0;
  n = 0;
  for(i=0; i<size[2]; i++)
    {
    for(j=0; j<size[1]; j++)
      {
      for(k=0; k<size[0]; k++)
        {
        ImageType::IndexType index;
        index[0]=k;
        index[1]=j;
        index[2]=i;
        if(!m_AEUseMask || m_AEMask->GetPixel(index)>0)
          {
          l = (m_Image->GetPixel(index)-m_RidgeOp->GetDataMin())
              / (m_RidgeOp->GetDataMax()-m_RidgeOp->GetDataMin())*255;
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
    m_AEThresh = (l+0.5)/256.0*(m_RidgeOp->GetDataMax()-m_RidgeOp->GetDataMin())+m_RidgeOp->GetDataMin();
    }
  else
    {
    m_AEThresh = m_RidgeOp->GetDataMax();
    }
  std::cout << "AEThesh = " << m_AEThresh << std::endl;
*/
}

/**
 * AutoExtract tubes */
template<class TInputImage, class TInputMask>   
bool
TubeNetExtractor<TInputImage, TInputMask>   
::AutoExtract(int , int ) // zMin, int zMax)
{
/*  itk::Size<3> size;
  size = m_Image->GetLargestPossibleRegion().GetSize();

  int i, j, k;
  for(i=zMin; i<zMax; i++)
    {
    for(j=5; j<size[1]-5; j++)
      {
      for(k=5; k<size[0]-5; k++)
        { 
        Index<3> index;
        index[0]=k;
        index[1]=j;
        index[2]=i;
        if( (!m_AEUseMask || m_AEMask->GetPixel(index)>0)
            && m_Image->GetPixel(index) > m_AEThresh
            && m_RidgeOp->GetDataMask()->GetPixel(index) == 0 )
          {
          ExtractTube(k, j, i);
          }
        }  

      if(j/10 == j/10.0 && m_StatusCallBack) 
        {
        char s[80];
        sprintf(s, "%d.%d of %d.%d", i, j, size[2]-1,size[1]-1);
        m_StatusCallBack("Auto Extract", s, 0);
        if(m_IdleCallBack)
          { 
          m_IdleCallBack();
          }
        }
      }
    }   
  if(m_StatusCallBack)
    {
    m_StatusCallBack("Auto Extract", "Done!", 1);
    }
*/ 
  return true;
}

/**
 * Generate vessel mask */
template<class TInputImage, class TInputMask>   
void
TubeNetExtractor<TInputImage, TInputMask>   
::DrawVesselMask(MaskType * ) // mask)
{
 /* // Needs to iterate
  TubeType::ChildrenListType * tubeList = m_TubeNet->GetChildren(0, "Tube");
  TubeType::ChildrenListType::iterator tubeIt;

  tubeIt = tubeList->begin();
  while(tubeIt != tubeList->end())
    {
    m_RidgeOp->DrawTube<MaskType>(mask, dynamic_cast<TubeType *>((*tubeIt).GetPointer()));
    ++tubeIt;
    }
  tubeList->clear();
*/
}

}; // end namespace itk

#endif
