/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkImageToTubeRigidMetric.h,v $
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
#ifndef __itkImageToTubeRigidMetric_h
#define __itkImageToTubeRigidMetric_h

#include "itkPoint.h"
#include "itkGroupSpatialObject.h"
#include "itkTubeSpatialObject.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkLinearInterpolateImageFunction.h"
#include <time.h>
#include "itkEuler3DTransform.h"
#include "itkImageToSpatialObjectMetric.h"
#include <itkGaussianDerivativeImageFunction.h>
//#include <itkGaussianSecondDerivativeImageFunction.h>

namespace itk
{
  
/** \class ImageToTubeRigidMetric
 * \brief Computes similarity between two objects to be registered
 */

template < class TFixedImage, class TMovingSpatialObject> 
class ITK_EXPORT ImageToTubeRigidMetric 
                           : public ImageToSpatialObjectMetric<TFixedImage,TMovingSpatialObject>
{
public:
  /** Standard "Self" typedef. */
  typedef ImageToTubeRigidMetric  Self;
  /** Standard "Superclass" typedef. */
  typedef ImageToSpatialObjectMetric<TFixedImage,TMovingSpatialObject>  Superclass;
  /** Smart pointer typedef support   */
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Type definition for a tube point */
  typedef TubeSpatialObjectPoint<3>  TubePointType;
  typedef TubeSpatialObject<3>       TubeType;
  typedef TMovingSpatialObject       MovingSpatialObjectType;
  typedef typename MovingSpatialObjectType::ChildrenListType ChildrenListType;
  typedef GroupSpatialObject<3>      TubeNetType;
  typedef Image<unsigned char,3>     MaskImageType;
  typedef typename MaskImageType::Pointer MaskImagePointer;
  typedef typename MaskImageType::IndexType   IndexType;
  typedef TFixedImage FixedImageType;
  typedef GaussianDerivativeImageFunction<TFixedImage> DerivativeImageFunctionType;
//  typedef GaussianSecondDerivativeImageFunction<TFixedImage> SecondDerivativeImageFunctionType;
  typedef typename Superclass::DerivativeType DerivativeType;
  typedef typename Superclass::ParametersType ParametersType;
  typedef typename Superclass::MeasureType MeasureType;

  typedef vnl_vector<double> VectorType;
  typedef vnl_matrix<double> MatrixType;
  typedef Point<double,3>    PointType;

  /** Run-time type information (and related methods).*/
  itkTypeMacro(ImageToTubeRigidMetric, ImageToSpatialObjectMetric);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Space dimension is the dimension of parameters space */
  enum { SpaceDimension = 6 }; //TMapper::SpaceDimension };
  enum { ImageDimension = 3 };
  enum { RangeDimension = 6 };

  unsigned int GetNumberOfParameters(void) const  {return SpaceDimension;};
  
  /** Typedef for the Range calculator */
  typedef MinimumMaximumImageCalculator<FixedImageType>   RangeCalculatorType;

  /** Type used for representing point components  */
  typedef typename Superclass::CoordinateRepresentationType 
                                                 CoordinateRepresentationType;

  /** Type definition for the size */
  typedef typename TFixedImage::SizeType   SizeType;

  /** Type definition for the pixel type */
  typedef typename TFixedImage::PixelType  PixelType;

  /**  Type of the Transform Base class */
  typedef Euler3DTransform<double> TransformType;

  typedef typename TransformType::Pointer            TransformPointer;
  typedef typename TransformType::InputPointType     InputPointType;
  typedef typename TransformType::OutputPointType    OutputPointType;
  typedef typename TransformType::ParametersType     TransformParametersType;
  typedef typename TransformType::JacobianType       TransformJacobianType;

  /** Get the Derivatives of the Match Measure */
  const DerivativeType & GetDerivative(const ParametersType & parameters) const;
  void GetDerivative( const ParametersType & parameters,
                                     DerivativeType & derivative ) const;

  /** Get the Value for SingleValue Optimizers */
  MeasureType    GetValue( const ParametersType & parameters ) const ;

  /** Get Value and Derivatives for MultipleValuedOptimizers */
  void GetValueAndDerivative( const ParametersType & parameters,
                    MeasureType & Value, DerivativeType  & Derivative ) const ;

  /** SubSample the MovingSpatialObject tube */
  void SubSampleTube(unsigned int sampling);

  /** Apply the center of rotation to the transformation*/
  ParametersType ApplyCenterOfRotation( const ParametersType & parameters );

  /** Set kappa value */
  itkSetMacro(Kappa,double);

  vnl_vector_fixed<double,3> GetCenterOfRotation(void) {return mC;}
 
  /** Apply a sparse registration to get closer */
  void SparseRegistration(ParametersType & parameters);

  /** Initialize the metric */
  void Initialize(void);

  /** Set the extent of the blurring */
  itkSetMacro(Extent,double);
  itkGetMacro(Extent,double);

  TransformPointer GetTransform(void) const {
    return dynamic_cast<TransformType*>(this->m_Transform.GetPointer());}

  itkSetObjectMacro(MaskImage, MaskImageType);

  itkSetObjectMacro(Verbose, bool);
  itkGetObjectMacro(Verbose, bool);
  
  itkSetObjectMacro(Sampling, unsigned int);
  itkGetObjectMacro(Sampling, unsigned int);

protected:

  ImageToTubeRigidMetric();
  virtual ~ImageToTubeRigidMetric();
  ImageToTubeRigidMetric(const Self&) {}
  void operator=(const Self&) {}

  void ComputeImageRange(void);
  void GetDeltaAngles(const Point<double,3> &  x,const vnl_vector_fixed<double,3> & dx, double *dA, double *dB, double *dG) const;

private:

  typename DerivativeImageFunctionType::Pointer m_DerivativeImageFunction;
 // typename SecondDerivativeImageFunctionType::Pointer m_SecondDerivativeImageFunction;
  MaskImagePointer                       m_MaskImage;
  typename ChildrenListType::iterator         TubeIterator;
  typename std::vector<TubePointType>::iterator  TubePointIterator;
  unsigned int                           m_NumberOfPoints;
  std::list<double>                      m_Weight;
  typename std::list<double>::iterator   WeightIterator;
  double                                 m_SumWeight;
  vnl_matrix<double>                     m_BiasV;
  vnl_matrix<double>                     m_BiasVI;
  double                                 m_ImageMin;
  double                                 m_ImageMax;
  typename RangeCalculatorType::Pointer  m_RangeCalculator;
  unsigned int                           m_Iteration;
  double                                 m_Kappa;
  double                                 m_RegImageThreshold;
  //MovingSpatialObjectPointer          m_MovingSpatialObject;
  //FixedImagePointer                   m_FixedImage;
  double                              m_Extent;
  mutable double                      m_Scale;
  unsigned int                        m_Sampling;

  mutable OutputPointType             m_CurrentPoint;
  mutable double                      m_BlurredValue;
  bool                                m_Verbose;
//  mutable VectorType                  m_Derivatives;

  double mScale;
  double mAlpha, mBeta, mGamma;
  long t1clock;

  vnl_matrix<double> * mT;
  //vnl_matrix<double> * mTI;
  vnl_vector<double> * mO;
  vnl_vector_fixed<double,3> mC;

  vnl_vector<double> * mTempV;
  vnl_vector<double> * mTempV2;

  //void SetTransform(vnl_matrix<double> * newT) const;
  void SetOffset(double oX, double oY, double oZ) const;
  void SetAngles(double newAlpha, double newBeta, double newGamma) const;
  void TransformPoint(vnl_vector<double> * in, vnl_vector<double> * out) const;
  void TransformVector(vnl_vector<double> * in, vnl_vector<double> * out);
  void TransformCoVector(vnl_vector<double> * in, vnl_vector<double> * out) const;

 
  /** Set the scale of the blurring */
  void SetScale(const double scale) const {m_Scale = scale;}
  itkGetConstMacro(Scale,double);

  /** Test whether the specified point is inside
   * Thsi method overload the one in the ImageMapper class
   * \warning This method cannot be safely used in more than one thread at
   * a time.
   * \sa Evaluate(); */
  bool IsInside( const InputPointType & point ) const;

  VectorType *  EvaluateAllDerivatives(void) const;
  VectorType *  ComputeThirdDerivatives(void) const;

  itkGetConstMacro(BlurredValue,double);

  VectorType* GetSecondDerivatives() const;
  MatrixType* GetHessian(PointType,double,double) const;
  double ComputeLaplacianMagnitude(Vector<double,3> *v) const;
  double  ComputeThirdDerivatives(Vector<double,3> *v) const;
  double  ComputeDerivatives(Vector<double,3> *v) const;

  bool IsInsideMask( const IndexType & index ) const;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToTubeRigidMetric.txx"
#endif

#endif


