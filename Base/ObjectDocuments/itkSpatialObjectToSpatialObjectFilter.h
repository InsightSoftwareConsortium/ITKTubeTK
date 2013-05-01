/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkSpatialObjectToSpatialObjectFilter.h,v $
  Language:  C++
  Date:      $Date: 2004/10/13 15:45:11 $
  Version:   $Revision: 1.4 $
  Author:    Julien Jomier

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkSpatialObjectToSpatialObjectFilter_h
#define __itkSpatialObjectToSpatialObjectFilter_h

#include "itkProcessObject.h"
#include "itkConceptChecking.h"

namespace itk
{

namespace tube
{

/** \class SpatialObjectToSpatialObjectFilter
 * \brief Base class for filters that take an image as input and produce an image as output.
 */
template <class TInputSpatialObject, class TOutputSpatialObject>
class ITK_EXPORT SpatialObjectToSpatialObjectFilter :
public ProcessObject
{
public:

  /** Standard class typedefs. */
  typedef SpatialObjectToSpatialObjectFilter<TInputSpatialObject,
                                             TOutputSpatialObject>  Self;
  typedef ProcessObject                                             Superclass;
  typedef SmartPointer<Self>                                        Pointer;
  typedef SmartPointer<const Self>                                  ConstPointer;

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
  const InputSpatialObjectType *  GetInput( void );
  const InputSpatialObjectType *  GetInput(unsigned int idx);

  void GenerateOutputInformation( void ) {} // do nothing
  void GenerateData( void ) {} // do nothing

protected:

  SpatialObjectToSpatialObjectFilter( void );
  ~SpatialObjectToSpatialObjectFilter( void );

  virtual void PrintSelf(std::ostream& os, Indent indent) const;

private:

  SpatialObjectToSpatialObjectFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

}; // End class SpatialObjectToSpatialObjectFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSpatialObjectToSpatialObjectFilter.txx"
#endif

#endif // End !defined(__itkSpatialObjectToSpatialObjectFilter_h)
