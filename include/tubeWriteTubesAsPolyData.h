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
#ifndef __tubeWriteTubesAsPolyData_h
#define __tubeWriteTubesAsPolyData_h

// ITK Includes
#include "itkProcessObject.h"
#include "itkGroupSpatialObject.h"
#include "itkTubeSpatialObject.h"

// TubeTK Includes
#include "tubeWrappingMacros.h"

namespace tube
{
/** \class WriteTubesAsPolyData
 *
 *  \ingroup TubeTK
 */

class WriteTubesAsPolyData:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef WriteTubesAsPolyData                         Self;
  typedef itk::ProcessObject                           Superclass;
  typedef itk::SmartPointer< Self >                    Pointer;
  typedef itk::SmartPointer< const Self >              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro( WriteTubesAsPolyData);

  /***/
  /***/
  /***/
  typedef itk::GroupSpatialObject<3>       GroupSpatialObjectType;
  typedef itk::TubeSpatialObject<3>        TubeSpatialObjectType;

  /** Set the source image. */
  void SetInput( GroupSpatialObjectType * grp )
  { m_GroupSpatialObject = grp; };

  void SetTubeInput( TubeSpatialObjectType * tube )
  { m_GroupSpatialObject = GroupSpatialObjectType::New();
    m_GroupSpatialObject->AddChild( tube ); };

  /** Set the filename for saving. */
  itkSetMacro( FileName, std::string );
  itkGetMacro( FileName, std::string );

  itkSetMacro( CenterlineFileName, std::string );
  itkGetMacro( CenterlineFileName, std::string );

  itkSetMacro( NumberOfSides, int );
  itkGetMacro( NumberOfSides, int );

  void Write( void )
  { this->Update(); };


protected:
  WriteTubesAsPolyData( void );
  ~WriteTubesAsPolyData() {};

  void Update() override;

  void PrintSelf( std::ostream & os, itk::Indent indent ) const override;

private:
  /** itktubeTubeExtractor parameters **/
  WriteTubesAsPolyData( const Self & );
  void operator=( const Self & );

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void SetInput( const DataObjectIdentifierType &, itk::DataObject * ) 
    override {};

  GroupSpatialObjectType::Pointer  m_GroupSpatialObject;

  std::string                      m_FileName;
  std::string                      m_CenterlineFileName;

  int                              m_NumberOfSides;
};

} // End namespace tube

#endif // End !defined( __tubeWriteTubesAsPolyData_h )
