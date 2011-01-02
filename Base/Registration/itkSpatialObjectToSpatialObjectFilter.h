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
#ifndef __itkSpatialObjectToSpatialObjectFilter_h
#define __itkSpatialObjectToSpatialObjectFilter_h

#include "itkProcessObject.h"
#include "itkConceptChecking.h"

namespace itk
{

/** \class SpatialObjectToSpatialObjectFilter
 * \brief Base class for filters that take an image as input and
 * produce an image as output.
 */
template <class TInputSpatialObject, class TOutputSpatialObject>
class ITK_EXPORT SpatialObjectToSpatialObjectFilter :
public ProcessObject
{
public:

  /** Standard class typedefs. */
  typedef SpatialObjectToSpatialObjectFilter<TInputSpatialObject,
    TOutputSpatialObject>                                  Self;
  typedef ProcessObject                                    Superclass;
  typedef SmartPointer<Self>                               Pointer;
  typedef SmartPointer<const Self>                         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SpatialObjectToSpatialObjectFilter,ProcessObject);

  /** Some convenient typedefs. */
  typedef TInputSpatialObject    InputSpatialObjectType;
  typedef TOutputSpatialObject   OutputSpatialObjectType;
  typedef typename InputSpatialObjectType::Pointer
                                 InputSpatialObjectPointer;
  typedef typename InputSpatialObjectType::ConstPointer
                                 InputSpatialObjectConstPointer;

  /** ImageDimension constants */
  itkStaticConstMacro(InputSpatialObjectDimension, unsigned int,
                      InputSpatialObjectType::ObjectDimension);
  itkStaticConstMacro(OutputSpatialObjectDimension, unsigned int,
                      OutputSpatialObjectType::ObjectDimension);

  /** Set/Get the image input of this process object.  */
  virtual void SetInput( const InputSpatialObjectType *object);
  virtual void SetInput( unsigned int, const InputSpatialObjectType * object);
  const InputSpatialObjectType *  GetInput(void);
  const InputSpatialObjectType *  GetInput(unsigned int idx);

  void GenerateOutputInformation() {}; // do nothing
  void GenerateData() {}; // do nothing

protected:
  SpatialObjectToSpatialObjectFilter();
  ~SpatialObjectToSpatialObjectFilter();

  virtual void PrintSelf(std::ostream& os, Indent indent) const;

private:
  SpatialObjectToSpatialObjectFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSpatialObjectToSpatialObjectFilter.txx"
#endif

#endif
