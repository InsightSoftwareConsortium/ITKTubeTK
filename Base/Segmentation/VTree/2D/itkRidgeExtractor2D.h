/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkRidgeExtractor2D.h,v $
  Language:  C++
  Date:      $Date: 2005/07/20 13:43:20 $
  Version:   $Revision: 1.10 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRidgeExtractor2D_h
#define __itkRidgeExtractor2D_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <cmath>
#include <list>

#include "itkSplineApproximation1D.h"
#include "itkSplineND.h"
#include "itkOptBrent1D.h"

#include "itkRadiusExtractor2D.h"
#include "itkBlurImageFunction.h"
#include <itkVesselTubeSpatialObject.h>


namespace itk 
{

/**
 * This class extract the ridge of a tube given an image
 * 
 * /sa itkRidgeExtractor2D
 * /todo -Implement the optimizer into itk
 *       -Use blur at a point using image function as soon as Josh give
 *        me the code.
 */

template <class TInputImage>             
class ITK_EXPORT RidgeExtractor2D : public Object 
{
public:

  /** 
   * Standard self typedef */
  typedef RidgeExtractor2D Self;
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
   * Type definition for the optimizer */
  typedef OptBrent1D OptimizerType;
  /*
   * Pointer type for the image */
  typedef typename TInputImage::Pointer  ImagePointer;
  /**
   * Const Pointer type for the image */
  typedef typename TInputImage::ConstPointer ImageConstPointer;
  /**
   * Type definition for the input image pixel type. */
  typedef typename TInputImage::PixelType PixelType;
  /**
   * Type definition for the size type. */
  typedef typename ImageType::SizeType SizeType;
  /** 
   * Type definition for the region type. */
  typedef typename ImageType::RegionType RegionType;
  /**
   * Type definition for the index type. */
  typedef typename ImageType::IndexType  IndexType; 

    /** Tube SpatialObject typedefs*/
  typedef VesselTubeSpatialObject<2> TubeType;
  typedef typename TubeType::Pointer TubePointer;
  typedef typename TubeType::TubePointType TubePointType;
  typedef typename TubeType::PointType PointType;
  
  /** Defines the type of vectors used */
  typedef vnl_vector<double> VnlVectorType; 
  typedef itk::Vector<double,2> VectorType;
  
  /** Defines the type of vectors used */
  typedef vnl_vector<int> VnlIntVectorType; 
  typedef itk::Vector<int,2> IntVectorType;

  /** Defines the type of matrix used */
  typedef vnl_matrix<double> VnlMatrixType; 
  typedef itk::Matrix<double,2> MatrixType;

  /**
   * Standard enum for the number of dimension  */
  enum {ImageDimension = ImageType::ImageDimension};
  /**
   * Set the input image */
  void SetInputImage(ImagePointer inputImage);
  /**
   * Get the input image */
  itkGetConstObjectMacro(Image,ImageType);
  /**
   * Get the mask image */
  itkGetConstObjectMacro(DataMask,ImageType);
  /**
   * Set Mode MR */
  void SetModeMR(bool ModeMR); 
  /**
   * Set Mode Retina */
  void SetModeRetina(bool ModeRetina); 
  /**
   * Set Mode CT */
  void SetModeCT(bool modeCT);
  /*  
   * Set Data Minimum */
  void SetDataMin(double dataMin); 
  /**
   * Get Data Minimum */
  itkGetMacro(DataMin,double);
  /*  
   * Set Data Maximum */
  void SetDataMax(double dataMax); 
  /**
   * Get Data Maximum */
  itkGetMacro(DataMax,double);
  /*  
   * Set Traversal Step size */
  itkSetMacro(StepX,double); 
  /**
   * Get Traversal Step size */
  itkGetMacro(StepX,double);
  /*  
   * Set Tangent change threshold */
  itkSetMacro(ThreshT,double); 
  /**
   * Get Tangent change threshold */
  itkGetMacro(ThreshT,double);
  /*  
   * Set Spatial change threshold */
  itkSetMacro(ThreshX,double); 
  /**
   * Get Spatial change threshold */
  itkGetMacro(ThreshX,double);
  /*  
   * Set P2Q2 Threshold */
  itkSetMacro(ThreshP2Q2,double); 
  /**
   * Get P2Q2 Threshold */
  itkGetMacro(ThreshP2Q2,double);
  /*  
   * Set P2Q2 start Threshold */
  itkSetMacro(ThreshP2Q2Start,double); 
  /**
   * Get P2Q2 start Threshold */
  itkGetMacro(ThreshP2Q2Start,double);
  /*  
   * Set EV  Threshold */
  itkSetMacro(ThreshEV,double); 
  /**
   * Get EV  Threshold */
  itkGetMacro(ThreshEV,double);
  /*  
   * Set EV ratio Threshold */
  itkSetMacro(ThreshEVRatio,double); 
  /**
   * Get EV ratio Threshold */
  itkGetMacro(ThreshEVRatio,double);
  /*  
   * Set EV ratio start Threshold */
  itkSetMacro(ThreshEVRatioStart,double); 
  /**
   * Get EV ratio start Threshold */
  itkGetMacro(ThreshEVRatioStart,double);
  /**
   * Set Extract Bound Minimum */
  void SetExtractBoundMin(IntVectorType* extractBoundMin);
  /**
   * get Extract Bound Minimum */
  IntVectorType* GetExtractBoundMin(void);
  /**
   * Set Extract Bound Maximum */
  void SetExtractBoundMax(IntVectorType* extractBoundMax);
  /**
   * get Extract Bound Maximum */
  IntVectorType* GetExtractBoundMax(void);
  /**
   * Get the data spline */
  SplineND* GetDataSpline(void);
  /**
   * Get the data spline 1D*/
  Spline1D* GetDataSpline1D(void);
  /**
   * Get the data spline optimizer*/
  Optimizer1D* GetDataSplineOptimizer(void);
  /**
   * Set the scale */
  void SetScale(double scale);
  /**
   * Get the scale */
  double GetScale(void);
  /**
   * Set the extent */ 
  void SetExtent(double extent);
  /**
   * Get the extent */
  double GetExtent(void);
  /**
   * Set the autoScale */
  void SetAutoScale(bool autoScale);
  /**
   * Get the autoScale */
  itkGetMacro(AutoScale,bool);
  /**
   * Set the dynamic Scale */
  void SetDynamicScale(bool dynamicScale);
  /**
   * Get the dynamicScale */
  itkGetMacro(DynamicScale,bool);
  /**
   * Get the autoScaleUsed */
  itkGetMacro(AutoScaleUsed,double);
  /**
   * Set the Recovery Maximum */
  itkSetMacro(RecoveryMax,int);
  /**
   * Get the Recovery Maximum */
  itkGetMacro(RecoveryMax,int);
  /**
   * Delete a tube */
  bool DeleteTube(TubeType * tube);
  /**
   * Add a tube */
  bool AddTube(TubeType * tube);
  /**
   * Set the radius Extractor */
  void  SetRadiusExtractor(RadiusExtractor2D<TInputImage> * radiusExtractor);
  /**
   * Compute the intensity at the point x */
  double  Intensity(VnlIntVectorType * x);
  /**the ridgness at point x */
  double  Ridgeness(VectorType * x);
  /**
   * Compute the local Ridge */
  bool   LocalRidge(VectorType * x);
  /**
   * Traverse the ridge one way */
  TubePointer  TraverseOneWay(VectorType * newX, VectorType * newT, VectorType  * newN, int dir);
 // void  TraverseOneWay(VectorType * newX, VectorType * newT, VectorType  * newN, int dir);
  /**
   * Extract */
  TubePointer  Extract(VectorType * x, int tubeID);
  /**
   * Smooth a tube */
  void   SmoothTubeX(TubeType * tube, int h);
  /**
   * Set the idle callback */
  void   IdleCallBack(bool (*idleCallBack)());
  /**
   * Set the status callback */
  void   StatusCallBack(void (*statusCallBack)(char *, char *, int));
  /**
   * Set if we should extract a ridge or a valley */
  itkSetMacro(ExtractRidge,bool);
  /**
   * Set if we should extract a ridge or a valley */
  void SetExtractValley(bool extractvalley)
  {
    m_ExtractRidge = !extractvalley;
  }

protected:
  RidgeExtractor2D();
  virtual ~RidgeExtractor2D();
  RidgeExtractor2D(const Self&) {}
  void operator=(const Self&) {}

private:

  ImagePointer    m_Image; 
   
  bool m_Debug;
  bool VERBOSE;

  typename BlurImageFunction<ImageType>::Pointer m_DataOp;

  ImagePointer     m_DataMask;

  bool             m_ExtractRidge;
  bool             m_ModeMR;   
  bool             m_AutoScale;
  bool             m_DynamicScale;
  double           m_AutoScaleUsed;
  double           m_DynamicScaleUsed;
  RadiusExtractor2D<TInputImage> * m_RadiusExtractor;

  int              m_RecoveryMax;
   
  double           m_DataMin;
  double           m_DataMax;
  double           m_DataRange;
   
  double           m_StepX;
  double           m_ThreshT;
  double           m_ThreshX;
  double           m_ThreshEV;

  IntVectorType m_ExtractBoundMin;
  IntVectorType m_ExtractBoundMax;
   
  SplineApproximation1D      m_DataSpline1D;
  OptimizerType       m_DataSplineOpt;
  SplineND *       m_DataSpline;
  UserFunc<vnl_vector<int> *, double> * m_SplineValueFunc;
   
  double           m_ThreshP2Q2;
  double           m_ThreshP2Q2Start;
  double           m_ThreshEVRatio;
  double           m_ThreshEVRatioStart;

   
 
  VectorType* m_X;
  double                   m_XP, m_XQ, m_XR;
  double                   m_XVal;

  VectorType* m_XD;
  VnlMatrixType* m_XH;
  VectorType m_XHEVal;
  VnlMatrixType* m_XHEVect;
   
  TubePointer     m_Tube;
  int         m_TubeID;
  int         m_TubePointCount;
  
  bool   m_CalcTangents(TubeType *);
  bool  (*m_IdleCallBack)();
  void  (*m_StatusCallBack)(char *, char *, int);

  // List to support tubepoints since push_front doesn't exist for std::vector
  std::list<TubePointType> m_TubePointList;

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRidgeExtractor2D.txx"
#endif

#endif /* __itkRidgeExtractor2D_h */


