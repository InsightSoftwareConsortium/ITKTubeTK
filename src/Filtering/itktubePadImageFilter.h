/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/
#ifndef __itktubePadImageFilter_h
#define __itktubePadImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkConceptChecking.h"

namespace itk {

namespace tube {

/** \class PadImageFilter
 * \brief Pad two images with zeros to make them suitable for a convolution
 * in the frequency domain.
 *
 * The requires two input images. The first one is supposed to be the image
 * to convolve and the second one the kernel.
 *
 * The filter pad with zeros the first image by the size of the second image.
 * The second image is centered and padded with zeros to have the same size
 * than the first output image.
 * After the transform, both images have the same LargestPossibleRegion.
 *
 * The option PadToPowerOfTwo can be set to true to force the size of the
 * images be a power of two - if the size of the padded image is 274 without
 * this option, it would be increased to 512 when PadToPowerOfTwo is true.
 * This option is makes the images usable with vnl's implementation of FFT.
 * PadToPowerOfTwo is false by default.
 *
 * \author Gaetan Lehmann
 *
 * \sa FFTShiftImageFilter NormalizeToConstantImageFilter
 *   FFTRealToComplexConjugateImageFilter
 */
template<class TInputImage, class TOutputImage=TInputImage>
class ITK_EXPORT PadImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef PadImageFilter Self;

  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;

  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                              InputImageType;
  typedef TOutputImage                             OutputImageType;
  typedef typename InputImageType::Pointer         InputImagePointer;
  typedef typename InputImageType::ConstPointer    InputImageConstPointer;
  typedef typename InputImageType::PixelType       InputImagePixelType;
  typedef typename OutputImageType::Pointer        OutputImagePointer;
  typedef typename OutputImageType::ConstPointer   OutputImageConstPointer;
  typedef typename OutputImageType::PixelType      OutputImagePixelType;
  typedef typename InputImageType::RegionType      RegionType;
  typedef typename InputImageType::IndexType       IndexType;
  typedef typename InputImageType::SizeType        SizeType;

  /** ImageDimension constants */
  itkStaticConstMacro( InputImageDimension, unsigned int,
                      TInputImage::ImageDimension );
  itkStaticConstMacro( OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension );
  itkStaticConstMacro( ImageDimension, unsigned int,
                      TOutputImage::ImageDimension );

  /** Standard New method. */
  itkNewMacro( Self );

  /** Runtime information support. */
  itkTypeMacro( PadImageFilter, ImageToImageFilter );

  /**
   * Set/Get whether the images must be padded to a size equal to a
   * power of two.
   * This is required for vnl implementation of FFT, but not for FFTW.
   * The default is false.
   */
  void SetPadToPowerOfTwo( bool v )
    {
    if( v )
      {
      this->SetGreatestPrimeFactor( 2 );
      }
    else
      {
      this->SetGreatestPrimeFactor( 13 );
      }
    }
  bool GetPadToPowerOfTwo() const
    {
    return m_GreatestPrimeFactor == 2;
    }
  itkBooleanMacro( PadToPowerOfTwo );

  /**
   * Set/Get the greatest prime factor allowed on the size of the padded
   * image. The filter increase the size of the image to reach a size with
   * the greatest prime factor smaller or equal to the specified value.
   * The default value is 13, which is the greatest prime number for which
   * the FFT are precomputed in FFTW, and thus gives very good performance.
   * A greatest prime factor of 2 produce a size which is a power of 2,
   * and thus is suitable for vnl base fft filters.
   * A greatest prime factor of 1 or less - typically 0 - disable the
   * extra padding.
   */
  itkGetConstMacro( GreatestPrimeFactor, int );
  itkSetMacro( GreatestPrimeFactor, int );

  /**
   * Set/Get the padding method.
   */
  typedef enum { NO_PADDING=0, ZERO_FLUX_NEUMANN=1, ZERO=2, MIRROR=3,
    WRAP=4 } PadMethod;
  itkGetConstMacro( PadMethod, int );
  itkSetMacro( PadMethod, int );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( InputHasPixelTraitsCheck,
    ( Concept::HasPixelTraits<InputImagePixelType> ) );
  itkConceptMacro( InputHasNumericTraitsCheck,
    ( Concept::HasNumericTraits<InputImagePixelType> ) );
  /** End concept checking */
#endif


protected:
  PadImageFilter();
  ~PadImageFilter() {};

  void PrintSelf( std::ostream& os, Indent indent ) const override;

  void GenerateInputRequestedRegion() override;
  void GenerateOutputInformation() override;

  /** Single-threaded version of GenerateData.  This filter delegates
   * to other filters. */
  void GenerateData() override;


private:
  PadImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  int m_GreatestPrimeFactor;
  int m_PadMethod;

  bool isPrime( int n )
    {
    int last = ( int )std::sqrt( static_cast<float>( n ) );
    for( int x=2; x<=last; x++ )
      {
      if( n%x == 0 )
        {
        return false;
        }
      }
    return true;
    }

  int greatestPrimeFactor( int n )
    {
    int v = 2;
    while( v <= n )
      {
      if( n%v == 0 && isPrime( v ) )
        {
        n /= v;
        }
      else
        {
        v += 1;
        }
      }
    return v;
    }

}; // end of class

} // end namespace tube

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubePadImageFilter.hxx"
#endif

#endif
