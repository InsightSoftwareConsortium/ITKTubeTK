/*=========================================================================

Library:   TubeTK

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

#ifndef __itktubeSpatialObjectToImageRegistrationMethod_h
#define __itktubeSpatialObjectToImageRegistrationMethod_h

#include "itkCommand.h"
#include "itkSpatialObject.h"
#include "itkDataObjectDecorator.h"

namespace itk
{

namespace tube
{

/** \class SpatialObjectToImageRegistrationMethod base class for the registration
 * methods.
 *
 * This class has a separate hierarchy from the ImageRegistrationMethod
 * defined in ITK.  The purpose of this class is to provide the common
 * functionalities of a registration method in a context that is easy to
 * use from the Registration Helper class that provides an even
 * higher-level, user-friendly interface to a generic image registration
 * problem.
 *
 */
template <unsigned int ObjectDimension, class TImage>
class SpatialObjectToImageRegistrationMethod
  : public ProcessObject
{

public:

  typedef SpatialObjectToImageRegistrationMethod Self;
  typedef ProcessObject                  Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;

  itkTypeMacro( SpatialObjectToImageRegistrationMethod, ProcessObject );

  itkNewMacro( Self );

  //
  // Custom Typedefs
  //
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TImage::ImageDimension );

  typedef Transform<double, ObjectDimension,
          itkGetStaticConstMacro( ImageDimension )>  TransformType;

  typedef DataObjectDecorator<TransformType>         TransformOutputType;

  typedef typename DataObject::Pointer               DataObjectPointer;
  typedef Superclass::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;

  typedef SpatialObject<ObjectDimension>             SpatialObjectType;
  typedef TImage                                     ImageType;

  typedef typename SpatialObjectType::PointType      SpatialObjectPointType;
  typedef typename ImageType::PointType              PointType;

  typedef SpatialObject<itkGetStaticConstMacro( ImageDimension )>
                                                     ImageMaskObjectType;

  typedef SpatialObject< ObjectDimension >
                                                     SpatialObjectMaskObjectType;

  //
  // Custom Methods
  //
  itkSetMacro( RegistrationNumberOfWorkUnits, unsigned int );
  itkGetMacro( RegistrationNumberOfWorkUnits, unsigned int );

  itkSetObjectMacro( Observer, Command );
  itkGetModifiableObjectMacro( Observer, Command );

  void SetFixedImage( const ImageType * fixedImage );
  itkGetConstObjectMacro( FixedImage, ImageType );

  void SetMovingSpatialObject( const SpatialObjectType * movingSpatialObject );
  itkGetConstObjectMacro( MovingSpatialObject, SpatialObjectType );

  void SetFixedImageMaskObject( const ImageMaskObjectType * maskObject );
  itkGetConstObjectMacro( FixedImageMaskObject, ImageMaskObjectType );
  itkSetMacro( UseFixedImageMaskObject, bool );
  itkGetMacro( UseFixedImageMaskObject, bool );

  void SetMovingSpatialObjectMaskObject(
    const SpatialObjectMaskObjectType * maskObject );
  itkGetConstObjectMacro( MovingSpatialObjectMaskObject,
    SpatialObjectMaskObjectType );
  itkSetMacro( UseMovingSpatialObjectMaskObject, bool );
  itkGetMacro( UseMovingSpatialObjectMaskObject, bool );

  itkSetMacro( ReportProgress, bool );
  itkGetMacro( ReportProgress, bool );
  itkBooleanMacro( ReportProgress );

  /** Return the output of the registration process, which is a Transform */
  const TransformOutputType * GetOutput( void ) const;

protected:

  SpatialObjectToImageRegistrationMethod( void );
  virtual ~SpatialObjectToImageRegistrationMethod( void );

  virtual void    Initialize( void );

  /** Method that actually computes the registration. This method is
   * intended to be overloaded by derived classes. Those overload,
   * however, must invoke this method in the base class. */
  void GenerateData( void ) override;

  /** Provide derived classes with access to the Transform member
   * variable. */
  itkSetObjectMacro( Transform, TransformType );
  itkGetModifiableObjectMacro( Transform, TransformType );

  void PrintSelf( std::ostream & os, Indent indent ) const override;

  using Superclass::MakeOutput;
  virtual DataObjectPointer  MakeOutput( DataObjectPointerArraySizeType
    idx ) override;

  ModifiedTimeType GetMTime( void ) const override;

protected:

  typename TransformType::Pointer m_Transform;

private:

  // Purposely not implemented
  SpatialObjectToImageRegistrationMethod( const Self & );
  // Purposely not implemented
  void operator =( const Self & );

  unsigned int m_RegistrationNumberOfWorkUnits;

  Command::Pointer m_Observer;

  typename ImageType::ConstPointer              m_FixedImage;
  typename SpatialObjectType::ConstPointer      m_MovingSpatialObject;

  bool                                          m_UseFixedImageMaskObject;
  typename ImageMaskObjectType::ConstPointer    m_FixedImageMaskObject;

  bool m_UseMovingSpatialObjectMaskObject;
  typename SpatialObjectMaskObjectType::ConstPointer
       m_MovingSpatialObjectMaskObject;

  bool m_ReportProgress;

};

} // tube

} // itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSpatialObjectToImageRegistrationMethod.hxx"
#endif

#endif
