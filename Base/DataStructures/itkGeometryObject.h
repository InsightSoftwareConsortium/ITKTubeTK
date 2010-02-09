/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGeometryObject.h,v $
  Language:  C++
  Date:      $Date: 2009-04-07 14:34:16 $
  Version:   $Revision: 1.68 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
 
#ifndef __itkGeometryObject_h 
#define __itkGeometryObject_h 
 
// Disable warning for long symbol names in this file only
#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <itkDataObject.h>
#include <itkObjectFactory.h>

namespace itk  
{ 

namespace Geometry
{

/** 
 * \class GeometryObject
 * \brief Base class for objects that implement 
 *
 * To implement your own spatial object, you need to derive from the
 * following class, which requires the definition of just a few pure
 * virtual functions.  Examples of such functions are ValueAt(),
 * IsEvaluableAt(), and IsInside(), each of which has a meaning
 * specific to each particular object type.
 */ 

template< unsigned int TDimension = 3> 
class GeometryObject 
  :public DataObject
{ 

public: 
  typedef GeometryObject<TDimension> Self;
  typedef DataObject                Superclass; 
  
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer; 

  /** Dimension of the object.  This constant is used by functions that are
   * templated over GeometryObject type when they need compile time access 
   * to the dimension of the object. */
  itkStaticConstMacro(ObjectDimension, unsigned int, TDimension);

  /** Run-time type information (and related methods). */ 
  itkTypeMacro( GeometryObject, DataObject );

  /** This defines the transformation from the global coordinate frame.
   *  By setting this transform, the local transform is computed */
  /*virtual void SetObjectTransform( TransformType * transform );
  itkGetObjectMacro(ObjectTransform,TransformType);
  itkGetConstObjectMacro(ObjectTransform,TransformType); */

  /** Copy information from the specified data set.  This method is
   * part of the pipeline execution model. By default, a ProcessObject
   * will copy meta-data from the first input to all of its
   * outputs. See ProcessObject::GenerateOutputInformation().  Each
   * subclass of DataObject is responsible for being able to copy
   * whatever meta-data it needs from from another DataObject.
   * ImageBase has more meta-data than its DataObject.  Thus, it must
   * provide its own version of CopyInformation() in order to copy the
   * LargestPossibleRegion from the input parameter. */
  virtual void CopyInformation(const DataObject *data);

  /** Update the information for this DataObject so that it can be used
   * as an output of a ProcessObject.  This method is used the pipeline
   * mechanism to propagate information and initialize the meta data
   * associated with a DataObject. This method calls its source's
   * ProcessObject::UpdateOutputInformation() which determines modified
   * times, LargestPossibleRegions, and any extra meta data like spacing,
   * origin, etc. */
  virtual void UpdateOutputInformation();

  /** Specify that the object has been updated */
  virtual void Update(void);

protected: 
 
  /** Constructor. */ 
  GeometryObject(); 

  /** Destructor. */ 
  virtual ~GeometryObject(); 

  virtual void PrintSelf( std::ostream& os, Indent indent ) const; 

private:

  GeometryObject(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
}; 

} // end of geometry namespace

} // end of namespace itk

#if !defined(CABLE_CONFIGURATION) 
#ifndef ITK_MANUAL_INSTANTIATION 
#include "itkGeometryObject.txx" 
#endif 
#endif 
 
#endif // __itkGeometryObject_h
