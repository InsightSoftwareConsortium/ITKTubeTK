/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __itkImageToImageDiffusiveDeformableRegistrationFilter_h
#define __itkImageToImageDiffusiveDeformableRegistrationFilter_h

#include "itkPDEDeformableRegistrationFilter.h"
#include "itkImageToImageDiffusiveDeformableRegistrationFunction.h"

namespace itk
{

/** \class itkImageToImageDiffusiveDeformableRegistrationFilter
 * \brief TODO
 *
 * TODO Insert more description + warnings here
 *
 * This class is templated over the fixed image type, moving image type and the
 * deformation field type.
 *
 * \sa itkImageToImageDiffusiveDeformableRegistrationFunction
 * \ingroup DeformableImageRegistration MultiThreaded
 */

template < class TFixedImage, class TMovingImage, class TDeformationField >
class ITK_EXPORT ImageToImageDiffusiveDeformableRegistrationFilter
  : public PDEDeformableRegistrationFilter< TFixedImage,
                                           TMovingImage,
                                           TDeformationField >
{
public:
  /** Standard class typedefs. */
  typedef ImageToImageDiffusiveDeformableRegistrationFilter   Self;
  typedef PDEDeformableRegistrationFilter< TFixedImage,
                                          TMovingImage,
                                          TDeformationField > Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Self, PDEDeformableRegistrationFilter);

protected:
  ImageToImageDiffusiveDeformableRegistrationFilter();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  // Purposely not implemented
  ImageToImageDiffusiveDeformableRegistrationFilter(const Self&);
  void operator=(const Self&); // Purposely not implemented
 
};

} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/itkImageToImageDiffusiveDeformableRegistrationFilter+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "itkImageToImageDiffusiveDeformableRegistrationFilter.txx"
#endif

#endif
