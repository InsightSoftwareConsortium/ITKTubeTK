/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkRidgeExtractor.h,v $
  Language:  C++
  Date:      $Date: 2005/06/02 20:30:55 $
  Version:   $Revision: 1.11 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRidgeExtractor_h
#define __itkRidgeExtractor_h

#include <cmath>
#include <list>

#include "itkSplineApproximation1D.h"
#include "itkSplineND.h"
#include "itkOptParabolicFit1D.h"
#include "itkRadiusExtractor.h"
#include "itkBlur3DImageFunction.h"

#include "itkContinuousIndex.h"
#include "itkVesselTubeSpatialObject.h"


namespace itk 
{

/**
 * This class extract the ridge of a tube given an image
 * 
 * /sa itkRidgeExtractor
 * /todo -Implement the optimizer into itk
 *       -Use blur at a point using image function as soon as Josh give
 *        me the code.
 */

template <class TInputImage>             
class ITK_EXPORT RidgeExtractor : public Object 
{
public:

  typedef RidgeExtractor             Self;
  typedef Object                     Superclass;
  typedef SmartPointer<Self>         Pointer;
  typedef SmartPointer<const Self>   ConstPointer;

  itkTypeMacro(RidgeExtractor, Object);

  itkNewMacro(RidgeExtractor);

  /**
   * Type definition for the input image. */
  typedef TInputImage                        ImageType;

  /**
   * Standard for the number of dimension */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      ::itk::GetImageDimension< TInputImage>::ImageDimension);

  /**
   * Type definition for the input image. */
  typedef Image<float, 3>                    MaskType;

  /**
   * Pointer type for the image */
  typedef typename TInputImage::Pointer      ImagePointer;

  /**
   * Const Pointer type for the image */
  typedef typename TInputImage::ConstPointer ImageConstPointer;

  /** Type definition for the input image pixel type. */
  typedef typename TInputImage::PixelType    PixelType;

  /** Type definition for the input image index type. */
  typedef typename TInputImage::IndexType    IndexType;

  /** Type definition for the input image index type. */
  typedef ContinuousIndex<double, 3> ContinuousIndexType;

  /** Defines the type of vectors used */
  typedef Vector<double, 3>                  VectorType;
  
  /** Defines the type of matrix used */
  typedef Matrix<double, 3, 3>               MatrixType;

  typedef Matrix<double, 3, 2>               NormalPlaneMatrixType;
  
  /** Tube SpatialObject typedefs*/
  typedef VesselTubeSpatialObject<3>         TubeType;
  typedef typename TubeType::Pointer         TubePointer;
  typedef typename TubeType::TubePointType   TubePointType;

  /**
   * Set the input image */
  void SetInputImage(ImagePointer inputImage);

  /**
   * Get the input image */
  itkGetObjectMacro(Image, ImageType);

  /**
   * Get the mask image */
  itkGetObjectMacro(DataMask, MaskType);

  /*  
   * Set Data Minimum */
  void SetDataMin(double dataMin); 

  /**
   * Get Data Minimum */
  itkGetMacro(DataMin, double);

  /**
   * Set Data Maximum */
  void SetDataMax(double dataMax); 

  /**
   * Get Data Maximum */
  itkGetMacro(DataMax, double);

  /**
   * Set Traversal Step size */
  itkSetMacro(StepX, double); 

  /**
   * Get Traversal Step size */
  itkGetMacro(StepX, double);

  /**
   * Set Tangent change threshold */
  itkSetMacro(ThreshT, double); 

  /**
   * Get Tangent change threshold */
  itkGetMacro(ThreshT, double);

  /**
   * Set Spatial change threshold */
  itkSetMacro(ThreshX, double); 

  /**
   * Get Spatial change threshold */
  itkGetMacro(ThreshX, double);

  /**
   * Set P2Q2 Threshold */
  itkSetMacro(ThreshP2Q2, double); 

  /**
   * Get P2Q2 Threshold */
  itkGetMacro(ThreshP2Q2, double);

  /**
   * Set P2Q2 start Threshold */
  itkSetMacro(ThreshP2Q2Start, double); 

  /**
   * Get P2Q2 start Threshold */
  itkGetMacro(ThreshP2Q2Start, double);

  /**
   * Set EV  Threshold */
  itkSetMacro(ThreshEV, double); 

  /**
   * Get EV  Threshold */
  itkGetMacro(ThreshEV, double);

  /**
   * Set EV ratio Threshold */
  itkSetMacro(ThreshEVRatio, double); 

  /**
   * Get EV ratio Threshold */
  itkGetMacro(ThreshEVRatio, double);

  /**
   * Set EV ratio start Threshold */
  itkSetMacro(ThreshEVRatioStart, double); 

  /**
   * Get EV ratio start Threshold */
  itkGetMacro(ThreshEVRatioStart, double);

  /**
   * Set Extract Bound Minimum */
  itkSetMacro(ExtractBoundMin, IndexType); 

  /**
   * get Extract Bound Minimum */
  itkGetMacro(ExtractBoundMin, IndexType); 

  /**
   * Set Extract Bound Maximum */
  itkSetMacro(ExtractBoundMax, IndexType); 

  /**
   * get Extract Bound Maximum */
  itkGetMacro(ExtractBoundMax, IndexType); 

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
   * Set the dynamic Scale */
  void SetDynamicScale(bool dynamicScale);

  /**
   * Get the dynamicScale */
  itkGetMacro(DynamicScale, bool);

  itkGetMacro(DynamicScaleUsed, double);

  /**
   * Set the Recovery Maximum */
  itkSetMacro(RecoveryMax, int);

  /**
   * Get the Recovery Maximum */
  itkGetMacro(RecoveryMax, int);

  /**
   * Delete a tube */
  bool DeleteTube(TubeType * tube);

  /**
   * Add a tube */
  bool AddTube(TubeType * tube);

  /**
   * Add a tube */
  template <class TDrawMask>
  void DrawTube(TDrawMask * drawMask, TubeType * tube);

  /**
   * Set the radius Extractor */
  void  SetRadiusExtractor(RadiusExtractor<TInputImage> * radiusExtractor);

  /**
   * Compute the intensity at the point x */
  double  Intensity(IndexType & x);

  /**the ridgness at point x */
  double  Ridgeness(ContinuousIndexType & x);
  /**
   * Compute the local Ridge */
  bool   LocalRidge(ContinuousIndexType & x);
  /**
   * Traverse the ridge one way */
  TubeType *  TraverseOneWay(ContinuousIndexType & newX, VectorType & newT,
                             NormalPlaneMatrixType & newN, int dir);
  /**
   * Extract */
  TubePointer  Extract(ContinuousIndexType & x, int tubeID);
  /**
   * Set the idle callback */
  void   IdleCallBack(bool (*idleCallBack)());
  /**
   * Set the status callback */
  void   StatusCallBack(void (*statusCallBack)(char *, char *, int));
 
  void SmoothTubeX(TubeType * tube, int h);

protected:

  RidgeExtractor();
  virtual ~RidgeExtractor();

  RidgeExtractor(const Self&) {}
  void operator=(const Self&) {}


private:

  ImagePointer    m_Image; 
   
  bool m_Debug;
  bool m_Verbose;

  typename Blur3DImageFunction<ImageType>::Pointer m_DataOp;

  MaskType::Pointer  m_DataMask;

  bool             m_DynamicScale;
  double           m_DynamicScaleUsed;
  RadiusExtractor<TInputImage> * m_RadiusExtractor;

  int              m_RecoveryMax;
   
  double           m_DataMin;
  double           m_DataMax;
  double           m_DataRange;
   
  double           m_StepX;
  double           m_ThreshT;
  double           m_ThreshX;
  double           m_ThreshEV;

  IndexType        m_ExtractBoundMin;
  IndexType        m_ExtractBoundMax;
   
  SplineApproximation1D                 m_DataSpline1D;
  OptParabolicFit1D                     m_DataSplineOpt;
  SplineND*                             m_DataSpline;
  UserFunc<vnl_vector<int> *, double> * m_SplineValueFunc;
   
  double           m_ThreshP2Q2;
  double           m_ThreshP2Q2Start;
  double           m_ThreshEVRatio;
  double           m_ThreshEVRatioStart;

  std::list<TubePointType>               m_TubePointList;
  
  ContinuousIndexType       m_X;
  double                    m_XP, m_XQ, m_XR;
  double                    m_XVal;

  VectorType       m_XD;
  MatrixType       m_XH;
  VectorType       m_XHEVal;
  MatrixType       m_XHEVect;
   
  TubePointer      m_Tube;
  int              m_TubeID;
  int              m_TubePointCount;
  
  bool  (*m_IdleCallBack)();
  void  (*m_StatusCallBack)(char *, char *, int);

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRidgeExtractor.txx"
#endif

#endif /* __itkRidgeExtractor_h */


