/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkRadiusExtractor.h,v $
  Language:  C++
  Date:      $Date: 2005/08/03 15:43:56 $
  Version:   $Revision: 1.12 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkRadiusExtractor_h
#define __itkRadiusExtractor_h

#include "itkOptParabolicFit1D.h"
#include "itkBlur3DImageFunction.h"
#include <itkVesselTubeSpatialObject.h>

#include <vxl_version.h>
#if VXL_VERSION_DATE_FULL > 20040406
# include <vnl/vnl_cross.h>
# define itk_cross_3d vnl_cross_3d
#else
# define itk_cross_3d cross_3d
#endif

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
class ITK_EXPORT RadiusExtractor : public Object 
{
public:

  /** 
   * Standard self typedef */
  typedef RadiusExtractor Self;
  typedef Object  Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  itkTypeMacro(RadiusExtractor, Object);
  itkNewMacro(RadiusExtractor);

  typedef VesselTubeSpatialObject<3> TubeType;
  typedef typename TubeType::Pointer TubePointer;
  typedef typename TubeType::TubePointType TubePointType;
  typedef typename TubeType::PointType PointType;

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
  typedef Vector<double, 3> VectorType; 

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
  OptParabolicFit1D & GetMedialnessOpt(void);

  /**
   * Compute Medialness and Branchness */    
  void     ComputeMnessBness(double pntR, double w,
                             double *kernPos, double *kernPosCnt,
                             double *kernNeg, double *kernNegCnt,
                             double *kernBrn, double *kernBrnCnt,
                             double &mness, double &bness,
                             bool doBNess);

  /**
   * Compute the medialness at a point */    
  double   MedialnessAtPoint(TubePointType pnt, double pntR,
                             bool doBNess=false,
                             bool newKern=true, 
                             double w=1);

  /**
   * Compute the medialness at a kernel */       
  double   MedialnessAtKern(std::list<TubePointType> * tube, double pntR, bool doBNess);

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

  RadiusExtractor();
  virtual ~RadiusExtractor();
  RadiusExtractor(const Self&) {}
  void operator=(const Self&) {}

private:

  ImagePointer    m_Image; 
  typename Blur3DImageFunction<ImageType>::Pointer m_DataOp;
  
  bool     m_Debug;
  bool     m_Verbose;
  
  double   m_DataMin;
  double   m_DataMax;
       
  int      m_NumRadiusPoints;
  int      m_RadiusPointSpacing;
       
  /** Determine if the algorithm extracts ridge or a valley */
  bool     m_ExtractRidge;

  double   m_Radius0;
  double   m_RadiusMin;
  double   m_RadiusMax;
  
  OptParabolicFit1D m_MedialnessOpt;
       
  double   m_ThreshWVal;
  double   m_ThreshWValStart;
      
  TubePointType                         * m_KernPntArray;
  std::vector<TubePointType>::iterator  * m_IterPntArray;
  int                                     m_ArrayLen;

  std::list<TubePointType>                m_Kern;
  

  double                     m_KernMedial;
  double                     m_KernBranch;
  UserFunc<double, double> * m_MedialnessAtKern;
      
  int                        m_KernNumT;
  double                     m_KernCosT[20], m_KernSinT[20];
  double                     m_KernPos[40], m_KernNeg[40];
  double                     m_KernPosCnt[40], m_KernNegCnt[40];
  double                     m_KernBrn[40];
  double                     m_KernBrnCnt[40];
  VectorType                 m_KernN0;
  VectorType                 m_KernN1;

  int                        m_TubePointCount;
  int                        m_TubeLength;

  double                     m_Scale;
  double                     m_Extent;   

  bool (*m_IdleCallBack)();
  void (*m_StatusCallBack)(char *, char *, int);

  void CalcKernArray(TubeType * tube);
  void CalcKernRadiiOneWay(int iStart, int iEnd, bool forward);
  void CalcKernMeasures(void);
  void ApplyKernMeasures(TubeType * tube);

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRadiusExtractor.txx"
#endif

#endif /* __itkRadiusExtractor_h */


