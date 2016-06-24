/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/
#ifndef __tubeShrinkWithBlendingImage_h
#define __tubeShrinkWithBlendingImage_h

#include "itktubeShrinkWithBlendingImageFilter.h"
#include "itkObject.h"
#include "tubeWrappingMacros.h"

namespace tube
{
/** \class ShrinkWithBlendingImage
 *
 *  \ingroup TubeTKITK
 */

template< typename TInputImage, typename TOutputImage >
class ShrinkWithBlendingImage:
  public itk::ProcessObject
{
public:
    /** Standard class typedefs. */
    typedef ShrinkWithBlendingImage                         Self;
    typedef itk::SmartPointer< Self >                       Pointer;
    typedef itk::SmartPointer< const Self >                 ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro( Self );

    /** Run-time type information (and related methods). */
    itkTypeMacro( ShrinkWithBlendingImage, Object );


    /** Typedef to images */
    typedef TOutputImage                          OutputImageType;
    typedef TInputImage                           InputImageType;
    typedef typename TInputImage::IndexType       InputIndexType;
    typedef typename TInputImage::SizeType        InputSizeType;

    itkStaticConstMacro( ImageDimension, unsigned int,
                         TInputImage::ImageDimension );

    itkStaticConstMacro( OutputImageDimension, unsigned int,
                         TOutputImage::ImageDimension );

    typedef itk::Vector< float, ImageDimension >       PointImagePixelType;
    typedef itk::Image< PointImagePixelType, OutputImageDimension >
                                                       PointImageType;

    typedef itk::FixedArray< unsigned int, ImageDimension > ShrinkFactorsType;
    typedef itk::tube::ShrinkWithBlendingImageFilter< InputImageType,
                            OutputImageType > FilterType;
    /** Set the shrink factors. Values are clamped to
     * a minimum value of 1. Default is 1 for all dimensions. */
    tubeWrapSetMacro(ShrinkFactors,ShrinkFactorsType,ShrinkWithBlendingFilter)
    void SetShrinkFactor(unsigned int i, unsigned int factor);
    unsigned int GetShrinkFactor(unsigned int i);

    tubeWrapSetMacro(NewSize,InputSizeType,ShrinkWithBlendingFilter);
    tubeWrapGetMacro(NewSize,InputSizeType,ShrinkWithBlendingFilter);

    /** Get/Set the shrink factors. */
    tubeWrapGetMacro(ShrinkFactors,ShrinkFactorsType,ShrinkWithBlendingFilter);

    tubeWrapSetMacro(Overlap,InputIndexType,ShrinkWithBlendingFilter);
    tubeWrapGetMacro(Overlap,InputIndexType,ShrinkWithBlendingFilter);

    tubeWrapSetMacro(BlendWithMean,bool,ShrinkWithBlendingFilter);
    tubeWrapGetMacro(BlendWithMean,bool,ShrinkWithBlendingFilter);

    tubeWrapSetMacro(BlendWithMax,bool,ShrinkWithBlendingFilter);
    tubeWrapGetMacro(BlendWithMax,bool,ShrinkWithBlendingFilter);

    tubeWrapSetMacro(BlendWithGaussianWeighting,bool,ShrinkWithBlendingFilter);
    tubeWrapGetMacro(BlendWithGaussianWeighting,bool,ShrinkWithBlendingFilter);

    tubeWrapSetMacro(UseLog,bool,ShrinkWithBlendingFilter);
    tubeWrapGetMacro(UseLog,bool,ShrinkWithBlendingFilter);

    tubeWrapGetObjectMacro(PointImage,PointImageType,ShrinkWithBlendingFilter);

    tubeWrapCallMacro(GenerateOutputInformation,ShrinkWithBlendingFilter);
    tubeWrapCallMacro(GenerateInputRequestedRegion,ShrinkWithBlendingFilter);


    tubeWrapSetConstObjectMacro(Input,InputImageType,ShrinkWithBlendingFilter);
    tubeWrapGetConstObjectMacro(Input,InputImageType,ShrinkWithBlendingFilter);
    tubeWrapCallMacro(Update,ShrinkWithBlendingFilter);

    tubeWrapGetObjectMacro(Output,OutputImageType,ShrinkWithBlendingFilter);

protected:
  ShrinkWithBlendingImage( void );
  ~ShrinkWithBlendingImage() {}
  void PrintSelf(std::ostream & os, itk::Indent indent) const;

private:
  /** itkShrinkWithBlendingImageFilter parameters **/
  ShrinkWithBlendingImage(const Self &);
  void operator=(const Self &);
  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput(const DataObjectIdentifierType &, itk::DataObject *) {};


  typename FilterType::Pointer m_ShrinkWithBlendingFilter;

};
} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeShrinkWithBlendingImage.hxx"
#endif

#endif // End !defined( __tubeShrinkWithBlendingImage_h )
