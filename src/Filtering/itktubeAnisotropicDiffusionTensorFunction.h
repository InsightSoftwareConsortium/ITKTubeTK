/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeAnisotropicDiffusionTensorFunction_h
#define __itktubeAnisotropicDiffusionTensorFunction_h

#include <itkDiffusionTensor3D.h>
#include <itkFiniteDifferenceFunction.h>
#include <itkImageRegionIterator.h>
#include <itkSymmetricSecondRankTensor.h>

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>

namespace itk
{

namespace tube
{

/** \class AnisotropicDiffusionTensorFunction
 * \brief This class is a function object that is used
 * to create a solver filter for edge enhancement diffusion equation
 *
 * \warning Does not handle image directions.  Re-orient images to axial
 * ( direction cosines = identity matrix ) before using this function.
 *
 * \sa AnisotropicDiffusionTensorImageFilter
 * \ingroup FiniteDifferenceFunctions
 * \ingroup Functions
 */
template< class TImageType >
class AnisotropicDiffusionTensorFunction
  : public FiniteDifferenceFunction< TImageType >
{
public:
  /** Standard class typedefs. */
  typedef AnisotropicDiffusionTensorFunction          Self;
  typedef FiniteDifferenceFunction<TImageType>        Superclass;
  typedef SmartPointer< Self >                        Pointer;
  typedef SmartPointer< const Self >                  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ) */
  itkTypeMacro( AnisotropicDiffusionTensorFunction,
    FiniteDifferenceFunction );

  /** Extract some parameters from the superclass. */
  itkStaticConstMacro( ImageDimension, unsigned int,
    Superclass::ImageDimension );

  /** Convenient typedefs. */
  typedef typename Superclass::TimeStepType             TimeStepType;
  typedef typename Superclass::PixelType                PixelType;
  typedef double                                        ScalarValueType;
  typedef typename Superclass::NeighborhoodType         NeighborhoodType;
  typedef typename Superclass::FloatOffsetType          FloatOffsetType;
  typedef typename Superclass::ImageType::SpacingType   SpacingType;

  /** Diffusion tensor typedefs. */
  typedef DiffusionTensor3D< double >              DiffusionTensorType;
  typedef itk::Image< DiffusionTensorType, 3 >     DiffusionTensorImageType;

  /** The default boundary condition for finite difference
   * functions that is used unless overridden in the Evaluate() method. */
  typedef ZeroFluxNeumannBoundaryCondition< DiffusionTensorImageType >
                                           DefaultBoundaryConditionType;
  typedef ConstNeighborhoodIterator< DiffusionTensorImageType,
    DefaultBoundaryConditionType>          DiffusionTensorNeighborhoodType;

  /** Scalar derivative typedefs. */
  typedef itk::Vector< ScalarValueType,
    itkGetStaticConstMacro( ImageDimension )>    ScalarDerivativeType;
  typedef itk::Image< ScalarDerivativeType, 3 >  ScalarDerivativeImageType;

  typedef ImageRegionIterator< ScalarDerivativeImageType >
    ScalarDerivativeImageRegionType;

  /** Tensor derivative typedefs. */
  typedef itk::Matrix< ScalarValueType,
    itkGetStaticConstMacro( ImageDimension ),
    itkGetStaticConstMacro( ImageDimension ) >  TensorDerivativeType;
  typedef itk::Image< TensorDerivativeType, 3 > TensorDerivativeImageType;
  typedef ImageRegionIterator< TensorDerivativeImageType >
    TensorDerivativeImageRegionType;

  /** A global data type for this class of equations.  Used to store
   * values that are needed in calculating the time step and other
   * intermediate products such as derivatives that may be used by
   * virtual functions called from ComputeUpdate.  Caching these values
   * here allows the ComputeUpdate function to be const and thread safe. */
  struct GlobalDataStruct
    {
    /** Hessian matrix
     * ( Second order partial derivatives of the intensity image ) */
    TensorDerivativeType  m_dxy;

    /** First order partial derivatives of the tensors */
    TensorDerivativeType  m_DT_dxy;

    /** First order partial derivatives of the intensity image */
    ScalarDerivativeType  m_dx;

    ScalarValueType       m_GradMagSqr;

    }; // End struct GlobalDataStruct

  /** Compute the equation value.  Inherited from the superclass: call
   *  one of the other two ComputeUpdate() functions instead. */
  virtual PixelType ComputeUpdate( const NeighborhoodType & neighborhood,
    void *globalData, const FloatOffsetType& = FloatOffsetType( 0.0 ) )
    override;

  /** Compute the equation value. The spacing of the images associated
   * with the given neighborhoods and regions should be the same as that
   * given. */
  virtual PixelType ComputeUpdate( const NeighborhoodType & neighborhood,
    const DiffusionTensorNeighborhoodType & tensorNeighborhood,
    const SpacingType & spacing, void * globalData,
    const FloatOffsetType & = FloatOffsetType( 0.0 ) );

  /** Compute the equation value using pre-computed derivatives. The
   * spacing of the images associated with the given neighborhoods and
   * regions should be the same as that given. */
  virtual PixelType ComputeUpdate(
    const DiffusionTensorNeighborhoodType & tensorNeighborhood,
    const ScalarDerivativeImageRegionType & intensityFirstDerivatives,
    const TensorDerivativeImageRegionType & intensitySecondDerivatives,
    const TensorDerivativeImageRegionType & tensorFirstDerivatives,
    void * globalData, const FloatOffsetType & = FloatOffsetType( 0.0 ) );

  /** Computes the time step for an update given a global data structure.
   *  Returns the time step supplied by the user. We don't need
   *  to use the global data supplied since we are returning a fixed value.
   *  */
  virtual TimeStepType ComputeGlobalTimeStep( void * itkNotUsed( globalData ) )
    const override
    { return this->GetTimeStep(); }

  /** Set/Get the time step. For this class of anisotropic diffusion
   * filters, the time-step is supplied by the user and remains fixed
   * for all updates. */
  void SetTimeStep( const TimeStepType & t )
    { m_TimeStep = t; }
  const TimeStepType & GetTimeStep( void ) const
    { return m_TimeStep; }

  /** Utility function to check whether the time step is stable, optionally
   * based on the spacing of the given image */
  template< class TPixel, unsigned int VImageDimension >
  void CheckTimeStepStability( const itk::Image< TPixel,
    VImageDimension > * input, bool useImageSpacing );

  /** Computes the first and second order partial derivatives of an
   * intensity image. */
  void ComputeIntensityFirstAndSecondOrderPartialDerivatives(
    const NeighborhoodType & neighborhood,
    ScalarDerivativeImageRegionType & firstOrderResult,
    TensorDerivativeImageRegionType & secondOrderResult,
    const SpacingType & spacing ) const;

  /** Computes the first order partial derivative of a diffusion tensor
   *  image. */
  void ComputeDiffusionTensorFirstOrderPartialDerivatives(
    const DiffusionTensorNeighborhoodType & tensorNeighborhood,
    TensorDerivativeImageRegionType & firstOrderResult,
    const SpacingType & spacing ) const;

  /** Determines whether to use the image spacing information in
   * calculations. Set the flag to ON if you want derivatives in physical
   * space, or OFF if you want derivatives in isotropic pixel space.
   * Default is ON. */
  void SetUseImageSpacing( bool newUseImageSpacing )
    { m_UseImageSpacing = newUseImageSpacing; }
  bool GetUseImageSpacing( void ) const
    { return m_UseImageSpacing; }

  /** Returns a pointer to a global data structure that is passed to this
   * object from the solver at each calculation. */
  virtual void * GetGlobalDataPointer( void ) const override
    {
    GlobalDataStruct *ans = new GlobalDataStruct();
    return ans;
    }

  virtual void ReleaseGlobalDataPointer( void * GlobalData ) const override
    { delete static_cast<GlobalDataStruct *>( GlobalData ); }

protected:
  AnisotropicDiffusionTensorFunction( void );
  virtual ~AnisotropicDiffusionTensorFunction( void ) {}

  void PrintSelf( std::ostream &s, Indent indent ) const override;

  /** The offset of the center pixel in the neighborhood. */
  unsigned int m_Center;

  /** Stride length along the y-dimension. */
  unsigned int m_xStride[itkGetStaticConstMacro( ImageDimension )];

  /** Defines various positions surrounding the center pixel in an image
    iterator. */
  unsigned int m_positionA[itkGetStaticConstMacro( ImageDimension )];
  unsigned int m_positionB[itkGetStaticConstMacro( ImageDimension )];
  unsigned int m_positionAa[itkGetStaticConstMacro( ImageDimension )]
      [itkGetStaticConstMacro( ImageDimension )];
  unsigned int m_positionBa[itkGetStaticConstMacro( ImageDimension )]
      [itkGetStaticConstMacro( ImageDimension )];
  unsigned int m_positionCa[itkGetStaticConstMacro( ImageDimension )]
      [itkGetStaticConstMacro( ImageDimension )];
  unsigned int m_positionDa[itkGetStaticConstMacro( ImageDimension )]
      [itkGetStaticConstMacro( ImageDimension )];

  /** Computes the first and second derivatives of an intensity image. */
  void ComputeIntensityFirstAndSecondOrderPartialDerivatives(
      const NeighborhoodType &neighborhood,
      ScalarDerivativeType &firstOrderResult,
      TensorDerivativeType &secondOrderResult,
      const SpacingType &spacing ) const;

  /** Computes the first order partial derivative of a diffusion tensor
   *  image. */
  void ComputeDiffusionTensorFirstOrderPartialDerivatives(
      const DiffusionTensorNeighborhoodType &tensorNeighborhood,
      TensorDerivativeType &firstOrderResult,
      const SpacingType &spacing ) const;

  /** Computes the final update term based on the results of the first and
    * second derivative computations */
  PixelType ComputeFinalUpdateTerm(
      const DiffusionTensorNeighborhoodType &tensorNeighborhood,
      const GlobalDataStruct* gd ) const;

private:
  //purposely not implemented
  AnisotropicDiffusionTensorFunction( const Self& );
  void operator=( const Self& );   //purposely not implemented

  TimeStepType    m_TimeStep;
  bool            m_UseImageSpacing;

}; // End class AnisotropicDiffusionTensorFunction

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeAnisotropicDiffusionTensorFunction.hxx"
#endif

#endif // End !defined( __itktubeAnisotropicDiffusionTensorFunction_h )
