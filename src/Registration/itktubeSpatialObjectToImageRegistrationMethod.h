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

#ifndef __itkSpatialObjectToImageRegistrationMethod_h
#define __itkSpatialObjectToImageRegistrationMethod_h

#include "itkSpatialObject.h"
#include "itkImageRegistrationMethod.h"

namespace itk
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
template <class TSpatialObject, class TImage>
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
  itkStaticConstMacro( SpatialObjectDimension, unsigned int,
                       TSpatialObject::Dimension );

  itkStaticConstMacro( ImageDimension, unsigned int,
                       TImage::ImageDimension );

  typedef Transform<double,
                    itkGetStaticConstMacro( SpatialObjectDimension ),
                    itkGetStaticConstMacro( ImageDimension )>
  TransformType;

  typedef DataObjectDecorator<TransformType> TransformOutputType;

  typedef typename DataObject::Pointer DataObjectPointer;
  typedef Superclass::DataObjectPointerArraySizeType
                                       DataObjectPointerArraySizeType;

  typedef TSpatialObject SpatialObjectType;
  typedef TImage         ImageType;

  typedef GroupSpatialObject<itkGetStaticConstMacro( SpatialObjectDimension )>
  GroupType;

  typedef typename TSpatialObject::PointType SpatialObjectPointType;
  typedef typename TImage::PointType PointType;

  typedef SpatialObject<itkGetStaticConstMacro( ImageDimension )>
  ImageMaskObjectType;

  typedef SpatialObject<itkGetStaticConstMacro( SpatialObjectDimension )>
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
  void SetMovingGroupSpatialObject( const GroupType * movingSpatialObject );

  itkGetConstObjectMacro( MovingGroupSpatialObject, GroupType );

  void SetRegionOfInterest( const PointType & point1,
    const PointType & point2 );

  itkSetMacro( UseRegionOfInterest, bool );
  itkGetMacro( UseRegionOfInterest, bool );
  itkSetMacro( RegionOfInterestPoint1, PointType );
  itkGetMacro( RegionOfInterestPoint1, PointType );
  itkSetMacro( RegionOfInterestPoint2, PointType );
  itkGetMacro( RegionOfInterestPoint2, PointType );

  void SetFixedImageMaskObject( const ImageMaskObjectType * maskObject );

  itkGetConstObjectMacro( FixedImageMaskObject, ImageMaskObjectType );

  itkSetMacro( UseFixedImageMaskObject, bool );
  itkGetMacro( UseFixedImageMaskObject, bool );

  void SetMovingSpatialObjectMaskObject(
    const SpatialObjectMaskObjectType * maskObject );

  itkGetConstObjectMacro( MovingImageMaskObject, SpatialObjectMaskObjectType );

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

  void PrintSelf( std::ostream & os, Indent indent ) const override;

  /** Provide derived classes with access to the Transform member
   * variable. */
  itkSetObjectMacro( Transform, TransformType );
  itkGetModifiableObjectMacro( Transform, TransformType );

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

  typename ImageType::ConstPointer       m_FixedImage;
  typename GroupType::ConstPointer       m_MovingGroupSpatialObject;

  bool      m_UseRegionOfInterest;
  PointType m_RegionOfInterestPoint1;
  PointType m_RegionOfInterestPoint2;

  bool                                   m_UseFixedImageMaskObject;
  typename MaskObjectType::ConstPointer  m_FixedImageMaskObject;

  bool                                   m_UseMovingSpatialObjectMaskObject;
  typename MaskObjectType::ConstPointer  m_MovingSpatialObjectMaskObject;

  bool m_ReportProgress;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSpatialObjectToImageRegistrationMethod.hxx"
#endif

#endif
