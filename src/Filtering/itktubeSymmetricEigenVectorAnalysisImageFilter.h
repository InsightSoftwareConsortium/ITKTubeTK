/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeSymmetricEigenVectorAnalysisImageFilter_h
#define __itktubeSymmetricEigenVectorAnalysisImageFilter_h

#include <itkSymmetricEigenAnalysis.h>
#include <itkUnaryFunctorImageFilter.h>

namespace itk
{

namespace tube
{

// This functor class invokes the computation of eigenanalysis for
// every pixel. The input pixel type must provide the API for the [][]
// operator, while the output pixel type must provide the API for the
// [] operator. Input pixel matrices should be symmetric.
//
// The default operation is to order eigenvalues in ascending order.
// You may also use OrderEigenValuesBy() to order eigenvalues by
// magnitude as is common with use of tensors in vessel extraction.
namespace Functor
{

#if !defined(ITK_LEGACY_REMOVE)
using OrderTypeOfEigenValue = itk::EigenValueOrderEnum;
#endif

template< class TInput, class TOutput, class TMatrix >
class SymmetricEigenVectorAnalysisFunction
{
public:
  SymmetricEigenVectorAnalysisFunction( void ) {}
  ~SymmetricEigenVectorAnalysisFunction( void ) {}
  typedef SymmetricEigenAnalysis< TInput, TOutput, TMatrix > CalculatorType;
  bool operator!=( const SymmetricEigenVectorAnalysisFunction & ) const
    {
    return false;
    }
  bool operator==( const SymmetricEigenVectorAnalysisFunction & other ) const
    {
    return !( *this != other );
    }

  inline TMatrix operator()( const TInput & x )
    {
    TOutput      eigenValues;
    eigenValues.Fill(0);
    TMatrix eigenVectorMatrix;
    m_Calculator.ComputeEigenValuesAndVectors( x, eigenValues,
      eigenVectorMatrix );
    return eigenVectorMatrix;
    }

  /** Method to explicitly set the dimension of the matrix */
  void SetDimension( unsigned int n )
    {
    m_Calculator.SetDimension( n );
    }

  unsigned int GetDimension() const
  {
    return m_Calculator.GetDimension();
  }

#if !defined(ITK_LEGACY_REMOVE)
  /** Enables reverse compatibility for enumeration values */
  using EigenValueOrderType = itk::EigenValueOrderEnum;
#endif
#if !defined(ITK_LEGACY_REMOVE)
  // We need to expose the enum values at the class level
  // for backwards compatibility
  static constexpr EigenValueOrderEnum OrderByValue = EigenValueOrderEnum::OrderByValue;
  static constexpr EigenValueOrderEnum OrderByMagnitude = EigenValueOrderEnum::OrderByMagnitude;
  static constexpr EigenValueOrderEnum DoNotOrder = EigenValueOrderEnum::DoNotOrder;
#endif

  /** Order eigenvalues. Default is to OrderByValue:  lambda_1 < lambda_2
   * < .... */
  void OrderEigenValuesBy( EigenValueOrderEnum order )
    {
    if( order == EigenValueOrderEnum::OrderByMagnitude )
      {
      m_Calculator.SetOrderEigenMagnitudes( true );
      }
    else if( order == EigenValueOrderEnum::DoNotOrder )
      {
      m_Calculator.SetOrderEigenValues( false );
      }
    }

private:
  CalculatorType m_Calculator;

}; // End class SymmetricEigenVectorAnalysisFunction

} // End namespace functor


/** \class SymmetricEigenVectorAnalysisImageFilter
 * \ingroup IntensityImageFilters  Multithreaded  TensorObjects
 */
template< class TInputImage, class TOutputImage, class TOutputMatrix >
class SymmetricEigenVectorAnalysisImageFilter
  : public UnaryFunctorImageFilter< TInputImage, TOutputMatrix,
    Functor::SymmetricEigenVectorAnalysisFunction<
      typename TInputImage::PixelType, typename TOutputImage::PixelType,
      typename TOutputMatrix::PixelType > >
{
public:
  /** Standard class typedefs. */
  typedef SymmetricEigenVectorAnalysisImageFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage, TOutputMatrix,
    Functor::SymmetricEigenVectorAnalysisFunction<
      typename TInputImage::PixelType, typename TOutputImage::PixelType,
      typename TOutputMatrix::PixelType > >        Superclass;

  typedef SmartPointer< Self >                     Pointer;
  typedef SmartPointer< const Self >               ConstPointer;

  typedef typename Superclass::OutputImageType     OutputImageType;
  typedef typename TOutputImage::PixelType         OutputPixelType;
  typedef typename TInputImage::PixelType          InputPixelType;
  typedef typename Superclass::FunctorType         FunctorType;

  /** Typedefs to order eigenvalues.
   * OrderByValue:      lambda_1 < lambda_2 < ....
   * OrderByMagnitude:  |lambda_1| < |lambda_2| < .....
   * DoNotOrder:        Default order of eigenvalues obtained after QL
   *                    method
   */
  //typedef typename FunctorType::EigenValueOrderEnum EigenValueOrderEnum;

  /** Order eigenvalues. Default is to OrderByValue:
   * lambda_1 < lambda_2 < .... */
  void OrderEigenValuesBy( EigenValueOrderEnum order )
    {
    this->GetFunctor().OrderEigenValuesBy( order );
    }

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Print internal ivars */
  void PrintSelf( std::ostream& os, Indent indent ) const override
    { this->Superclass::PrintSelf( os, indent ); }

  /** Set the dimension of the tensor. ( For example the
   * SymmetricSecondRankTensor is a pxp matrix ) */
  void SetDimension( unsigned int p )
    {
    this->GetFunctor().SetDimension( p );
    }

protected:
  SymmetricEigenVectorAnalysisImageFilter( void ) {}
  virtual ~SymmetricEigenVectorAnalysisImageFilter( void ) {}

private:
  //purposely not implemented
  SymmetricEigenVectorAnalysisImageFilter( const Self& );
  //purposely not implemented
  void operator=( const Self& );

}; // End class SymmetricEigenVectorAnalysisImageFilter

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeSymmetricEigenVectorAnalysisImageFilter_h )
