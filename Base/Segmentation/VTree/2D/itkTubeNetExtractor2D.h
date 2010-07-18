/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkTubeNetExtractor2D.h,v $
  Language:  C++
  Date:      $Date: 2003/01/13 19:59:27 $
  Version:   $Revision: 1.4 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTubeNetExtractor2D_h
#define __itkTubeNetExtractor2D_h

#include "itkTubeExtractor.h"
#include <itkTubeNet.h>

namespace itk 
{

/**
 * This class extract the a tube given an image
 * 
 * /sa itkTubeNetExtractor2D
 * /todo Implement the optimizer into itk
 *       Use blur at a point using image function as soon as Josh give
 *           me the code.
 */

template <class TInputImage>             
class ITK_EXPORT TubeNetExtractor2D : public TubeExtractor<TInputImage> 
{
public:

  /** 
   * Standard self typedef */
  typedef TubeNetExtractor2D Self;
  /**
   * Standard "Superclass" typedef. */
  typedef Object  Superclass;
  /** 
   * Smart pointer typedef support */
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  /** 
   * Run-time type information (and related methods). */
  itkTypeMacro(Self, Object);
  /**
   * Method for creation through the object factory. */
  itkNewMacro(Self);
  /**
   * Type definition for the input image. */
  typedef TInputImage  ImageType;
  /**
   * Pointer type for the image */
  typedef typename TInputImage::Pointer  ImagePointer;
  /**
   * Const Pointer type for the image */
  typedef typename TInputImage::ConstPointer ImageConstPointer;
  /**
   * Type definition for the input image pixel type. */
  typedef typename TInputImage::PixelType PixelType;
  /**
	 * Defines the type of vectors used
	 */
	typedef vnl_vector<double> VectorType; 
  /**
	 * Standard enum for the number of dimension
	 */
  enum {ImageDimension = ImageType::ImageDimension};
  /**
   * Set the input image */
  void SetInputImage(ImagePointer inputImage);
  /**
   * Get the input image */
  itkGetConstObjectMacro(Image,ImageType);
  /**
   * Extract a 2D tube */
  bool ExtractTube(float x, float y);
  /**
   * Delete a tube */
  bool DeleteTube(Tube * newTube);
  /**
   * Get the tube Net */
  TubeNet::Pointer GetTubeNet(void);
  /**
   * Set use mask */
  itkSetMacro(AEUseMask,bool);
  /**
   * Get use mask */
  itkGetMacro(AEUseMask,bool);
  /**
   * Set the mask  */    
  void SetAutoExtractMask(ImagePointer autoExtractMask);
  /**
   * Get the mask  */   
  itkGetConstObjectMacro(AEMask,ImageType);
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
   * Save a tube net */
  //bool   Save(char *fname);
  /**
   * Load a tube net*/
  //bool   Load(char *fname);  
  /**
   * Set the tube callback */
  void   NewTubeCallBack(void (*newTubeCallBack)(Tube *));

protected:
  TubeNetExtractor2D();
  virtual ~TubeNetExtractor2D();
  TubeNetExtractor2D(const Self&) {}
  void operator=(const Self&) {}

private:

  ImagePointer    m_Image; 
   
  TubeNet::Pointer  m_TubeNet;
  int               m_TubeNum;

  bool                    m_AEUseMask;
  ImagePointer            m_AEMask;
  float                  *m_AEThresh;


};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeNetExtractor2D.txx"
#endif

#endif /* __itkTubeNetExtractor2D_h */


