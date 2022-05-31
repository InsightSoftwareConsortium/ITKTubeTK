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

#ifndef __itktubeSpatialObjectSource_h
#define __itktubeSpatialObjectSource_h

#include <itkProcessObject.h>
#include <itkSpatialObject.h>

namespace itk
{

namespace tube
{
/** \class SpatialObjectSource
 *
 * \brief Base class for all ProcessObject's that output SpatialObject's.
 *
 * SpatialObjectSource is the base class for all process objects that output
 * SpatialObject data. Specifically, this class defines the GetOutput() method
 * that returns a pointer to the output SpatialObject.
 */
template< class TOutputSpatialObject >
class SpatialObjectSource : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef SpatialObjectSource        Self;
  typedef ProcessObject              Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  typedef TOutputSpatialObject       OutputSpatialObjectType;

  typedef SpatialObject< TOutputSpatialObject::ObjectDimension >
    SpatialObjectType;
  typedef Superclass::DataObjectIdentifierType
    DataObjectIdentifierType;
  typedef Superclass::DataObjectPointerArraySizeType
    DataObjectPointerArraySizeType;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( SpatialObjectSource, ProcessObject );

  OutputSpatialObjectType * GetOutput( void );
  const OutputSpatialObjectType * GetOutput( void ) const;

  OutputSpatialObjectType * GetOutput( unsigned int idx );

  /** Graft the specified DataObject onto this ProcessObject's output.
   * This method grabs a handle to the specified DataObject's bulk
   * data to used as its output's own bulk data. It also copies the
   * region ivars ( RequestedRegion, BufferedRegion,
   * LargestPossibleRegion ) and meta-data ( Spacing, Origin ) from the
   * specified data object into this filter's output data object. Most
   * importantly, however, it leaves the Source ivar untouched so the
   * original pipeline routing is intact. This method is used when a
   * process object is implemented using a mini-pipeline which is
   * defined in its GenerateData() method.  The usage is:
   *
   * \code
   *    // setup the mini-pipeline to process the input to this filter
   *    firstFilterInMiniPipeline->SetInput( this->GetInput() );

   *    // setup the mini-pipeline to calculate the correct regions
   *    // and write to the appropriate bulk data block
   *    lastFilterInMiniPipeline->GraftOutput( this->GetOutput() );
   *
   *    // execute the mini-pipeline
   *    lastFilterInMiniPipeline->Update();
   *
   *    // graft the mini-pipeline output back onto this filter's output.
   *    // this is needed to get the appropriate regions passed back.
   *    this->GraftOutput( lastFilterInMiniPipeline->GetOutput() );
   * \endcode
   *
   * For proper pipeline execution, a filter using a mini-pipeline
   * must implement the GenerateInputRequestedRegion(),
   * GenerateOutputRequestedRegion(), GenerateOutputInformation() and
   * EnlargeOutputRequestedRegion() methods as necessary to reflect
   * how the mini-pipeline will execute ( in other words, the outer
   * filter's pipeline mechanism must be consistent with what the
   * mini-pipeline will do ).
   *  */
  virtual void GraftOutput( DataObject *output );

  /** Graft the specified data object onto this ProcessObject's named
   * output. This is similar to the GraftOutput method except it
   * allows you to specify which output is affected.
   * See the GraftOutput for general usage information.
   */
  virtual void GraftOutput( const DataObjectIdentifierType & key,
    DataObject *output );

  /** Graft the specified data object onto this ProcessObject's idx'th
   * output. This is similar to the GraftOutput method except it
   * allows you to specify which output is affected. The specified index
   * must be a valid output number ( less than
   * ProcessObject::GetNumberOfIndexedOutputs() ). See the GraftOutput for
   * general usage information. */
  virtual void GraftNthOutput( unsigned int idx, DataObject *output );
  using Superclass::MakeOutput;
  virtual ProcessObject::DataObjectPointer
    MakeOutput( ProcessObject::DataObjectPointerArraySizeType idx ) override;

protected:
  SpatialObjectSource( void );
  virtual ~SpatialObjectSource( void ) {}

private:
  SpatialObjectSource( const Self & ); // purposely not implemented
  void operator=( const Self & );      // purposely not implemented

}; // End class SpatialObjectSource

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeSpatialObjectSource.hxx"
#endif

#endif // End !defined( __itktubeSpatialObjectSource_h )
