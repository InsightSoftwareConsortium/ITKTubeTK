/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkTubeExtractor2D.h,v $
  Language:  C++
  Date:      $Date: 2005/07/20 13:43:20 $
  Version:   $Revision: 1.7 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkTubeExtractor2D_h
#define __itkTubeExtractor2D_h

#include "itkRidgeExtractor2D.h"
#include "itkRadiusExtractor2D.h"
#include "itkVesselTubeSpatialObject.h"


namespace itk 
{

/**
 * This class extract the a tube given an image
 * 
 * /sa itkTubeExtractor2D
 * /todo Implement the optimizer into itk
 *       Use blur at a point using image function as soon as Josh give
 *           me the code.
 */

template <class TInputImage>             
class ITK_EXPORT TubeExtractor2D : public Object 
{
public:

  /** 
   * Standard self typedef */
  typedef TubeExtractor2D Self;
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
  
  /**  Type definition for TubeSpatialObjec*/
  typedef VesselTubeSpatialObject<2> TubeType;
  typedef typename TubeType::Pointer TubePointer;
  typedef typename TubeType::TubePointType TubePointType;
  /**
   * Defines the type of vectors used
   */
  //typedef vnl_vector<double> VectorType; 
  typedef itk::Vector<double,2> VectorType;


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
  /*  
   * Set Data Minimum */
  void SetDataMin(double dataMin); 
  /**
   * Get Data Minimum */
  double GetDataMin(void); 
  /*  
   * Set Data Maximum */
  void SetDataMax(double dataMax); 
  /**
   * Get Data Maximum */
  double GetDataMax(void); 
  /*  
   * Set the radius */
  void SetRadius(double radius); 
  /**
   * Get the radius */
  double GetRadius(void); 
  /*  
   * Set Extract Ridge */
  void ExtractRidge(bool extractRidge); 
  /**
   * Get Extract Ridge */
  bool ExtractRidge(void); 
  /*  
   * Set Extract Valley */
  void ExtractValley(bool extractValley); 
  /**
   * Get Extract Valley */
  bool ExtractValley(void); 
  /**
   * Get the ridge extractor */
  typename RidgeExtractor2D<TInputImage>::Pointer GetRidgeOp(void);
  /**
   * Get the radius extractor */
  typename RadiusExtractor2D<TInputImage>::Pointer GetRadiusOp(void);
  /**
   * Extract the 2D tube given the position of the first point
   * and the tube ID */
  bool ExtractTube(float x, float y, unsigned int tubeID); 
  /**
   * Delete a tube */
  bool DeleteTube(TubeType* tube);
  /** 
   * Get the last tube extracted */
  TubePointer GetLastTube(void);
  /**
   * Get the last position */
  VectorType* GetLastPosition(void);
  /**
   * Set the tube color */
  void SetColor(float* color);
  /**
   * Get the tube color */
  itkGetMacro(Color,float*);
  /**
   * Set the idle callback */
  void   IdleCallBack(bool (*idleCallBack)());
  /**
   * Set the status callback */
  void   StatusCallBack(void (*statusCallBack)(char *, char *, int));
  /**
   * Set the tube callback */
  void   NewTubeCallBack(void (*newTubeCallBack)(Tube *));
  /**
   * Set the status callback */
  void   AbortProcess(bool (*abortProcess)());
  /**
   * Set the mode for RoP */
  void SetModeRetina(bool modeRetina);

protected:
  TubeExtractor2D();
  virtual ~TubeExtractor2D();
  TubeExtractor2D(const Self&) {}
  void operator=(const Self&) {}

  typename RidgeExtractor2D<ImageType>::Pointer  m_RidgeOp;
  typename RadiusExtractor2D<ImageType>::Pointer m_RadiusOp;
   TubePointer  m_Tube;
  bool (*m_IdleCallBack)();
  void (*m_StatusCallBack)(char *, char *, int);
  void (*m_NewTubeCallBack)(TubeType *);
  bool (*m_AbortProcess)();

private:

  bool m_Debug;
  ImagePointer    m_Image;   
  VectorType *m_X;
  float* m_Color;

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTubeExtractor2D.txx"
#endif

#endif /* __itkTubeExtractor2D_h */


