/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkRadiusExtractor2D.h,v $
  Language:  C++
  Date:      $Date: 2005/09/27 22:02:14 $
  Version:   $Revision: 1.9 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRadiusExtractor2D_h
#define __itkRadiusExtractor2D_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

//#include "itkOptParabolicFit1D.h"
#include "itkOptGoldenMean1D.h"
#include <itkTube.h>
#include "itkBlurImageFunction.h"
#include "itkVesselTubeSpatialObject.h"

namespace itk 
{

/**
 * This class extract the radius of a tube given an image
 * 
 * /sa itkRidgeExtractor
 * /todo Implement the optimizer into itk
 *       Use blur at a point using image function as soon as Josh give
 *           me the code.
 */

template <class TInputImage>             
class ITK_EXPORT RadiusExtractor2D : public Object 
{
public:

  /** 
   * Standard self typedef */
  typedef RadiusExtractor2D Self;

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
   * Type definition for the optimizer*/
  typedef OptGoldenMean1D OptimizerType;
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
   * Defines the type of vectors used */
  typedef vnl_vector<double> VectorType; 


  /** Tube spatial object typedefs */
  typedef VesselTubeSpatialObject<2> TubeType;
  typedef typename TubeType::Pointer TubePointer;
  typedef typename TubeType::TubePointType TubePointType;
  typedef typename TubeType::PointType PointType;
  typedef typename TubeType::VectorType ITKVectorType;

  /** 
   * Standard enum for the number of dimension */
  enum {ImageDimension = ImageType::ImageDimension};
  /** 
   * Set the input image */
  void SetInputImage(ImagePointer inputImage);
  /**
   * Get the input image */
  itkGetConstObjectMacro(Image,ImageType);
  /**
   * Set Scale */
  void SetScale(double scale);
  /**
   * Get Scale */
  itkGetMacro(Scale,double);
  /**
   * Set Extent */
  void SetExtent(double extent);
  /**
   * Get Extent */
  itkGetMacro(Extent,double);
  /**
   * Set Mode MR */
  itkSetMacro(ModeMR,bool); 
  /**
   * Set Mode CT */
  void SetModeCT(bool modeCT);
  /**
   * Set Mode Retina */
  void SetModeRetina(bool modeRetina);
  /*  
   * Set Data Minimum */
  itkSetMacro(DataMin,double); 
  /**
   * Get Data Minimum */
  itkGetMacro(DataMin,double);
  /*  
   * Set Data Maximum */
  itkSetMacro(DataMax,double); 
  /**
   * Get Data Maximum */
  itkGetMacro(DataMax,double);
  /**
   * Set Minimum Radius */
  void SetRadiusMin(double radiusMin);
  /**
   * Get Minimum Radius */
  itkGetMacro(RadiusMin,double);  
  /**
   * Set Maximum Radius */
  void SetRadiusMax(double radiusMax);
  /**
   * Get Maximum Radius */
  itkGetMacro(RadiusMax,double);
  /**
   * Set Radius0 */
  itkSetMacro(Radius0,double);
  /**
   * Get Radius0 */
  itkGetMacro(Radius0,double); 
  /**
   * Set ThreshWVal */
  itkSetMacro(ThreshWVal,double);
  /**
   * Get ThreshWVal */
  itkGetMacro(ThreshWVal,double); 
  /**
   * Set ThreshWVal Start */
  itkSetMacro(ThreshWValStart,double);
  /**
   * Get ThreshWVal Start*/
  itkGetMacro(ThreshWValStart,double); 
  /**
   * Set Extract Ridge */
  itkSetMacro(ExtractRidge,bool);
  /**
   * Get ExtractRidge*/
  itkGetMacro(ExtractRidge,bool); 
  /**
   * Return the optimizer */
  OptimizerType & GetMedialnessOpt(void);

  /**
   * Compute Medialness and Branchness */    
 /* void     ComputeMnessBness(double pntR, double w,
                             double *kernPos, double *kernPosCnt,
                             double *kernNeg, double *kernNegCnt,
                             double *kernBrn, double *kernBrnCnt,
                             double *mness, double *bness,
                             bool doBNess);*/
  /**
   * Compute the medialness at a point */    
  /*double   MedialnessAtPoint(TubePoint * pnt, double pntR,
                             bool doBNess=false,
                             bool newKern=true, 
                             double w=1);*/

  double MedialnessAtPoint(TubePointType pnt,double pntR);

  double BranchnessAtPoint(TubePointType pnt,double pntR);
 
  /**
   * Compute the medialness at a kernel */
  double   MedialnessAtKern(std::list<TubePointType> * tube, double pntR);
  //double   MedialnessAtKern(Tube * tube, double pntR, bool doBNess);

  /**
   * Compute the Branchness at a kernel */
  double  BranchnessAtKern(TubeType * tube, double pntR);

  /**
   * Calculate the optimal scale */    
  bool  CalcOptimalScale(TubePointType pnt, bool firstGuess =false);
  /**
   * Calculate Radii one way */    
  bool  CalcRadiiOneWay(std::vector<TubePointType>::iterator tubePntFrom,
                       std::vector<TubePointType>::iterator tubePntTo,
                      bool forward=true);
  
  /**
   * Calculate Radii */    
  bool     CalcRadii(TubeType * tube);
       
  void     SetIdleCallBack(bool (*idleCallBack)());
  void     SetStatusCallBack(void (*statusCallBack)(char *, char *, int));

protected:
  RadiusExtractor2D();
  virtual ~RadiusExtractor2D();
  RadiusExtractor2D(const Self&) {}
  void operator=(const Self&) {}

private:

  ImagePointer    m_Image; 
  typename BlurImageFunction<ImageType>::Pointer m_DataOp;
  
  bool m_Debug;

  
  double   m_DataMin;
  double   m_DataMax;
       
  bool     m_ModeMR;
       
  int      m_NumRadiusPoints;
  int      m_RadiusPointSpacing;
       
  /** Determine if the algorithm extracts ridge or a valley */
  bool     m_ExtractRidge;

  double   m_Radius0;
  double   m_RadiusMin;
  double   m_RadiusMax;
  
  OptimizerType m_MedialnessOpt;
       
  double   m_ThreshWVal;
  double   m_ThreshWValStart;
      
  TubePoint ** m_KernPntArray;
  std::list<TubePoint *>::iterator * m_IterPntArray;
  int      m_ArrayLen;
  std::list<TubePointType>*    m_Kern;
  double   m_KernMedial;
  double   m_KernBranch;
  UserFunc<double, double> * m_MedialnessAtKern;
      
  int m_KernNumT;
  double m_KernCosT[20], m_KernSinT[20];
  double m_KernPos[40], m_KernNeg[40];
  double m_KernPosCnt[40], m_KernNegCnt[40];
  double m_KernBrn[40];
  double m_KernBrnCnt[40];
  VectorType* m_KernN0;
  VectorType* m_KernN1;

  int m_TubePointCount;
  int m_TubeLength;

  double m_Scale;
  double m_Extent;   
  bool (*m_IdleCallBack)();
  void (*m_StatusCallBack)(char *, char *, int);

  //void CalcKernArray(Tube * tube);
  //void CalcKernRadiiOneWay(int iStart, int iEnd, bool forward);
  //void CalcKernMeasures(void);
  //void ApplyKernMeasures(Tube * tube);

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRadiusExtractor2D.txx"
#endif

#endif /* __itkRadiusExtractor2D_h */


