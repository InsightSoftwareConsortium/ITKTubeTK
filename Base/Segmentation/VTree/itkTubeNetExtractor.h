/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkTubeNetExtractor.h,v $
  Language:  C++
  Date:      $Date: 2003/12/16 14:07:06 $
  Version:   $Revision: 1.5 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTubeNetExtractor_h
#define __itkTubeNetExtractor_h

#include "itkTubeExtractor.h"
//#include <TubeIO/itkTubeNet.h>

namespace itk 
{

/**
 * This class extract the a tube given an image
 * 
 * /sa itkTubeNetExtractor
 * /todo Implement the optimizer into itk
 *       Use blur at a point using image function as soon as Josh gives
 *           me the code.
 */

template <class TInputImage, class TInputMask>             
class ITK_EXPORT TubeNetExtractor : public TubeExtractor<TInputImage> 
{
public:

  /** 
   * Standard self typedef */
  typedef TubeNetExtractor          Self;
  typedef Object                    Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  itkTypeMacro(Self, Object);

  itkNewMacro(Self);

  /**
   * Type definition for the input image. */
  typedef TInputImage  ImageType;

  /**
   * Type definition for the input mask. */
  typedef TInputMask                       MaskType;

  /**
   * Pointer type for the image */
  typedef typename ImageType::Pointer      ImagePointer;

  /**
   * Pointer type for the mask */
  typedef typename MaskType::Pointer       MaskPointer;

  /**
   * Input image pixel type. */
  typedef typename TInputImage::PixelType  PixelType;

  /**
   * Input image index type. */
  typedef typename TInputImage::IndexType  IndexType;

  /** Typedef for TubeSpatialObject */
  typedef VesselTubeSpatialObject<3>             TubeType;
  typedef typename TubeType::Pointer       TubePointer;

  /**
   * Defines the type of vectors used
   */
  typedef Vector<double, 3>                VectorType; 

  /**
   * Standard for the number of dimension
   */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      ::itk::GetImageDimension< TInputImage>::ImageDimension);

  /**
   * Set the input image */
  void SetInputImage(ImagePointer inputImage);

  /**
   * Get the input image */
  itkGetConstObjectMacro(Image, ImageType);

  /**
   * Extract a 3D tube */
  bool ExtractTube(float x, float y, float z);

  /**
   * Delete a tube */
  bool DeleteTube(TubeType * newTube);

  /**
   * Get the tube Net */
  TubeType::Pointer GetTubeNet(void);

  /**
   * Set use mask */
  itkSetMacro(AEUseMask,bool);

  /**
   * Get use mask */
  itkGetMacro(AEUseMask,bool);

  /**
   * Set the mask  */    
  void SetAutoExtractMask(MaskPointer autoExtractMask);

  /**
   * Get the mask  */   
  itkGetConstObjectMacro(AEMask, MaskType);

  /**
   * Auto extract tubes using a mask */    
  double  AutoExtractThresh(void);

  /**
   * Auto extract tubes using a mask */      
  void   AutoExtractThresh(double newAEThresh);

  /**
   * Auto extract tubes using a mask */      
  void   AutoExtractAutoThresh(double alpha=0.002);

  /**
   * AutoExtract tubes */
  bool   AutoExtract(int zMin, int zMax);

  /**
   * Generate output mask */
  void   DrawVesselMask(MaskType * maskImage);

  /**
   * Save a tube net */
  //bool   Save(char *fname);

  /**
   * Load a tube net*/
  //bool   Load(char *fname);  

  /**
   * Set the tube callback */
  void   NewTubeCallBack(void (*newTubeCallBack)(TubeType *));

protected:
  TubeNetExtractor();
  virtual ~TubeNetExtractor();
  TubeNetExtractor(const Self&) {}
  void operator=(const Self&) {}

private:

  ImagePointer      m_Image; 
  TubePointer       m_TubeNet;
  int               m_TubeNum;

  bool              m_AEUseMask;
  MaskPointer       m_AEMask;
  float             m_AEThresh;

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeNetExtractor.txx"
#endif

#endif /* __itkTubeNetExtractor_h */
