/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ITKHeader.h,v $
  Language:  C++
  Date:      $Date: 2007-07-10 11:35:36 -0400 ( Tue, 10 Jul 2007 ) $
  Version:   $Revision: 0 $

  Copyright ( c ) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkImageToImageRegistrationMethod_h
#define __itkImageToImageRegistrationMethod_h

#include "itkSpatialObject.h"
#include "itkImageRegistrationMethod.h"

namespace itk
{

/** \class ImageToImageRegistrationMethod base class for the registration
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
template <class TImage>
class ImageToImageRegistrationMethod
  : public ProcessObject
{

public:

  typedef ImageToImageRegistrationMethod Self;
  typedef ProcessObject                  Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;

  itkTypeMacro( ImageToImageRegistrationMethod, ProcessObject );

  itkNewMacro( Self );

  //
  // Custom Typedefs
  //
  itkStaticConstMacro( ImageDimension, unsigned int,
                       TImage::ImageDimension );

  typedef Transform<double,
                    itkGetStaticConstMacro( ImageDimension ),
                    itkGetStaticConstMacro( ImageDimension )>
  TransformType;

  typedef DataObjectDecorator<TransformType> TransformOutputType;

  typedef typename DataObject::Pointer DataObjectPointer;
  typedef Superclass::DataObjectPointerArraySizeType
                                       DataObjectPointerArraySizeType;

  typedef TImage ImageType;

  typedef typename TImage::PointType PointType;

  typedef SpatialObject<itkGetStaticConstMacro( ImageDimension )>
  MaskObjectType;

  //
  // Custom Methods
  //
  itkSetMacro( RegistrationNumberOfWorkUnits, unsigned int );
  itkGetMacro( RegistrationNumberOfWorkUnits, unsigned int );

  itkSetObjectMacro( Observer, Command );
  itkGetObjectMacro( Observer, Command );

  void SetFixedImage( const ImageType * fixedImage );

  itkGetConstObjectMacro( FixedImage, ImageType );

  void SetMovingImage( const ImageType * movingImage );

  itkGetConstObjectMacro( MovingImage, ImageType );

  void SetRegionOfInterest( const PointType & point1,
    const PointType & point2 );

  itkSetMacro( UseRegionOfInterest, bool );
  itkGetMacro( UseRegionOfInterest, bool );
  itkSetMacro( RegionOfInterestPoint1, PointType );
  itkGetMacro( RegionOfInterestPoint1, PointType );
  itkSetMacro( RegionOfInterestPoint2, PointType );
  itkGetMacro( RegionOfInterestPoint2, PointType );

  void SetFixedImageMaskObject( const MaskObjectType * maskObject );

  itkGetConstObjectMacro( FixedImageMaskObject, MaskObjectType );

  itkSetMacro( UseFixedImageMaskObject, bool );
  itkGetMacro( UseFixedImageMaskObject, bool );

  void SetMovingImageMaskObject( const MaskObjectType * maskObject );

  itkGetConstObjectMacro( MovingImageMaskObject, MaskObjectType );

  itkSetMacro( UseMovingImageMaskObject, bool );
  itkGetMacro( UseMovingImageMaskObject, bool );

  itkSetMacro( ReportProgress, bool );
  itkGetMacro( ReportProgress, bool );
  itkBooleanMacro( ReportProgress );

  /** Return the output of the registration process, which is a Transform */
  const TransformOutputType * GetOutput( void ) const;

protected:

  ImageToImageRegistrationMethod( void );
  virtual ~ImageToImageRegistrationMethod( void );

  virtual void    Initialize( void );

  /** Method that actually computes the registration. This method is
   * intended to be overloaded by derived classes. Those overload,
   * however, must invoke this method in the base class. */
  void GenerateData( void ) override;

  void PrintSelf( std::ostream & os, Indent indent ) const override;

  /** Provide derived classes with access to the Transform member
   * variable. */
  itkSetObjectMacro( Transform, TransformType );
  itkGetObjectMacro( Transform, TransformType );
  itkGetConstObjectMacro( Transform, TransformType );

  using Superclass::MakeOutput;
  virtual DataObjectPointer  MakeOutput( DataObjectPointerArraySizeType
    idx ) override;

  ModifiedTimeType GetMTime( void ) const override;

protected:

  typename TransformType::Pointer m_Transform;

private:

  // Purposely not implemented
  ImageToImageRegistrationMethod( const Self & );
  // Purposely not implemented
  void operator =( const Self & );

  unsigned int m_RegistrationNumberOfWorkUnits;

  Command::Pointer m_Observer;

  typename ImageType::ConstPointer       m_FixedImage;
  typename ImageType::ConstPointer       m_MovingImage;

  bool      m_UseRegionOfInterest;
  PointType m_RegionOfInterestPoint1;
  PointType m_RegionOfInterestPoint2;

  bool                                   m_UseFixedImageMaskObject;
  typename MaskObjectType::ConstPointer  m_FixedImageMaskObject;

  bool                                   m_UseMovingImageMaskObject;
  typename MaskObjectType::ConstPointer  m_MovingImageMaskObject;

  bool m_ReportProgress;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkImageToImageRegistrationMethod.hxx"
#endif

#endif
