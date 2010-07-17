/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkImageToTubeRigidRegistration.h,v $
  Language:  C++
  Date:      $Date: 2006/06/12 18:34:27 $
  Version:   $Revision: 1.1 $
  Author:    Julien Jomier

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkImageToTubeRigidRegistration_h
#define __itkImageToTubeRigidRegistration_h

#include "itkImageToSpatialObjectRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageToTubeRigidMetric.h"
#include "itkGradientDescentVariableStepOptimizer.h"
#include "itkGradientDescentOptimizer.h"
#include "itkImage.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "itkConjugateGradientOptimizer.h"
#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include "itkEuler3DTransform.h"
#include "itkImageRegionIterator.h"
#include "itkVectorContainer.h"
#include "itkTubeSpatialObject.h"
#include "itkCommandIterationUpdate.h"

namespace itk
{

/** \class ImageToTubeRigidRegistration
 * \brief Base class for registration methods
 *
 * This Class define the generic interface for a registration method.
 * The basic elements of a registration method are:
 *   - Metric to compare the FixedImage and the TMovingTube
 *   - Transformation used to register the FixedImage against the TMovingTube
 *   - Optimization method used to search for the best transformation
 * 
 * Registration is not limited to Images, and for this reason
 * this class is templated over the type of the FixedImage object,
 * the TMovingTube object and the transformation. This types are obtained
 * from the Metric type, to reduce the number of redundant
 * template parameters
 *
 *  \ingroup AffineImageRegistration 
 */

template <class TFixedImage, class TMovingTube>
class ITK_EXPORT ImageToTubeRigidRegistration 
: public ImageToSpatialObjectRegistrationMethod<TFixedImage,TMovingTube>
{
public:
  /** Standard "Self" typedef. */
  typedef ImageToTubeRigidRegistration  Self;

  /** Standard "Superclass" typedef. */
  typedef ImageToSpatialObjectRegistrationMethod<TFixedImage,TMovingTube> Superclass;
  typedef typename Superclass::FixedImageType FixedImageType;
  
  /** Smart pointer typedef support */
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Typedef of the mask image */
  typedef Image<unsigned char,3>          MaskImageType;
  typedef typename MaskImageType::Pointer MaskImagePointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToTubeRigidRegistration, 
               ImageToSpatialObjectRegistrationMethod);

  /**  Type of the metric. */
  typedef ImageToTubeRigidMetric<FixedImageType,TMovingTube>  MetricType;

  /**  Type of the parameters. */
  typedef typename MetricType::TransformParametersType    ParametersType;

  /**  Type of the transformation. */
  typedef typename MetricType::TransformType TransformType;
  typedef typename TransformType::Pointer    TransformPointer;
  /**  Dimension of the images.  */
  enum {ImageDimension = FixedImageType::ImageDimension,
    ParametersDimension = TransformType::ParametersDimension}; 

  /** Type of the Interpolator */
  typedef typename Superclass::InterpolatorType InterpolatorType;
  typedef typename InterpolatorType::Pointer InterpolatorPointer;

  /**  Type of the optimizer.  */
  //typedef GradientDescentVariableStepOptimizer  OptimizerType;
  typedef GradientDescentOptimizer  OptimizerType;
  //typedef OnePlusOneEvolutionaryOptimizer OptimizerType;
 
  /** Typedef for the optimizer observer */
  typedef itk::CommandIterationUpdate<  
                                  OptimizerType >    CommandIterationType;

  /** Get Center of Rotation */
  itk::Vector<double,3> GetCenterOfRotation(void)
  {
    itk::Vector<double,3> centerOfRotation;
    typename MetricType::Pointer metric = dynamic_cast<MetricType*>(this->GetMetric());
    for(unsigned int i=0;i<3;i++)
    {
      centerOfRotation[i]=metric->GetCenterOfRotation()(i);
    }
    return centerOfRotation;
  }

  /** Method that initiates the registration. */
  void StartRegistration(void);

  /** Start the sparse registration */
  void SparseRegistration(ParametersType & parameters);

  /** Set the number of iteration */
  itkSetMacro(NumberOfIteration,unsigned int);

  /** Set the learning rate */
  itkSetMacro(LearningRate,double);

  /** Set the initial position */
  void SetInitialPosition(double position[6]);

  /** Set the parameters scales */
  void SetParametersScale(double scales[6]);

  /** Set the iteration observer */
  itkSetObjectMacro(IterationCommand,CommandIterationType);

  /** Initialize the registration */
  void Initialize() throw (ExceptionObject);

  /** Set the mask image */
  itkSetObjectMacro(MaskImage,MaskImageType);

  /** Write a matlab file to display the metric */
  void CreateMatlabMetric(const char* filename);

  void SetExtent(float extent) {m_Extent = extent;}
  void SetVerbose(bool verbose) {m_Verbose = verbose;}
  void SetKappa(float kappa) {m_Kappa = kappa;}
  void SetSampling(unsigned int sampling) {m_Sampling = sampling;}

protected:

  ImageToTubeRigidRegistration();
  virtual ~ImageToTubeRigidRegistration() {};

private:

  ImageToTubeRigidRegistration(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MaskImagePointer              m_MaskImage;
  unsigned int                  m_NumberOfIteration;
  CommandIterationType::Pointer m_IterationCommand;
  bool                          m_IsInitialized;
  double                        m_LearningRate;
  ParametersType                m_InitialPosition;
  ParametersType                m_ParametersScale;
  float                         m_Extent;
  float                         m_Kappa;
  bool                          m_Verbose;
  unsigned int                  m_Sampling;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToTubeRigidRegistration.txx"
#endif

#endif
