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

#ifndef __itktubeTubeXIO_hxx
#define __itktubeTubeXIO_hxx

#include <cmath>

#include "itktubeTubeXIO.h"

#include "tubeStringUtilities.h"
#include "metaUtils.h"

namespace itk
{

namespace tube
{

template< unsigned int TDimension >
TubeXIO< TDimension >
::TubeXIO()
{
  m_TubeGroup = TubeGroupType::New();
  for( unsigned int i = 0; i < TDimension; ++i )
    {
    m_Dimensions[i] = 1;
    }
}

template< unsigned int TDimension >
TubeXIO< TDimension >
::~TubeXIO( void )
{
}

//
template< unsigned int TDimension >
void
TubeXIO< TDimension >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( this->m_TubeGroup.IsNotNull() )
    {
    os << indent << "Tube Group = " << this->m_TubeGroup << std::endl;
    }
  else
    {
    os << indent << "Tube Group = NULL" << std::endl;
    }
}

template< unsigned int TDimension >
bool
TubeXIO< TDimension >
::Read( const std::string & _fileName )
{
  std::vector< MET_FieldRecordType * > fields;

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "NDims", MET_INT, false );
  fields.push_back( mF );
  int nDimsRecNum = MET_GetFieldRecordNumber( "NDims", &fields );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "Dimensions", MET_FLOAT_ARRAY, false, nDimsRecNum );
  fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "VoxelSize", MET_FLOAT_ARRAY, false, nDimsRecNum );
  fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "NObjects", MET_INT, false );
  fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "End Header", MET_NONE, false );
  mF->terminateRead = true;
  fields.push_back( mF );

  std::ifstream tmpReadStream( _fileName.c_str(), std::ios::binary |
    std::ios::in );

  if( !tmpReadStream.rdbuf()->is_open() )
    {
    for( unsigned int i=0; i<fields.size(); ++i )
      {
      delete fields[i];
      }
    fields.clear();
    return false;
    }

  std::vector< MET_FieldRecordType * > extraFields;
  if( !MET_Read( tmpReadStream, &fields, ':', false, true, &extraFields ) )
    {
    std::cerr << "Tube: Read failed" << std::endl;
    for( unsigned int i=0; i<fields.size(); ++i )
      {
      delete fields[i];
      }
    fields.clear();
    return false;
    }

  mF = MET_GetFieldRecord( "NDims", &fields );
  if( mF->defined )
    {
    int nDims = ( int )mF->value[0];

    if( ( int )mF->value[0] != TDimension )
      {
      std::cerr << "Tube: Read failed: object is " << nDims
        << " dimensional and was expecting " << TDimension << " dimensional."
        << std::endl;
      for( unsigned int i=0; i<fields.size(); ++i )
        {
        delete fields[i];
        }
      fields.clear();
      return false;
      }
    }

  mF = MET_GetFieldRecord( "Dimensions", &fields );
  if( mF->defined )
    {
    for( unsigned int i=0; i<TDimension; ++i )
      {
      this->m_Dimensions[i] = ( int )( mF->value[i] );
      }
    }

  double voxelSize[ TDimension ];
  mF = MET_GetFieldRecord( "VoxelSize", &fields );
  if( mF->defined )
    {
    for( unsigned int i=0; i<TDimension; ++i )
      {
      voxelSize[i] = ( float )( mF->value[i] );
      }
    }
  else
    {
    for( unsigned int i=0; i<TDimension; ++i )
      {
      voxelSize[i] = 1;
      }
    }

  int nObjects = 0;
  mF = MET_GetFieldRecord( "NObjects", &fields );
  if( mF->defined )
    {
    nObjects = ( int )mF->value[0];
    }

  for( unsigned int i=0; i<fields.size(); ++i )
    {
    delete fields[i];
    }
  fields.clear();

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "ID", MET_INT, false );
  fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "Type", MET_STRING, false );
  fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "Anat", MET_STRING, false );
  fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "TreeType", MET_STRING, false );
  fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "VParent", MET_INT, false );
  fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "Attachpt", MET_INT, false );
  fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "Color", MET_STRING, false );
  fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "PointDim", MET_STRING, true );
  fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "NPoints", MET_INT, true );
  fields.push_back( mF );

  mF = new MET_FieldRecordType;
  MET_InitReadField( mF, "Points", MET_NONE, true );
  mF->terminateRead = true;
  fields.push_back( mF );

  for( int obj = 0; obj < nObjects; ++obj )
    {
    if( !MET_Read( tmpReadStream, &fields, ':', false, true, &extraFields ) )
      {
      std::cerr << "Tube: Read failed" << std::endl;
      for( unsigned int i=0; i<fields.size(); ++i )
        {
        delete fields[i];
        }
      fields.clear();
      return false;
      }

    int tubeId = 0;
    mF = MET_GetFieldRecord( "ID", &fields );
    if( mF->defined )
      {
      tubeId = ( int )( mF->value[0] );
      }

    char tubeType[80];
    strcpy( tubeType, "" );
    mF = MET_GetFieldRecord( "Type", &fields );
    if( mF->defined )
      {
      strcpy( tubeType, ( char * )( mF->value ) );
      }

    char tubeAnat[80];
    strcpy( tubeAnat, "" );
    mF = MET_GetFieldRecord( "Anat", &fields );
    if( mF->defined )
      {
      strcpy( tubeAnat, ( char * )( mF->value ) );
      }

    char tubeTreeType[80];
    strcpy( tubeTreeType, "" );
    mF = MET_GetFieldRecord( "TreeType", &fields );
    if( mF->defined )
      {
      strcpy( tubeTreeType, ( char * )( mF->value ) );
      }

    int tubeParentId = -1;
    mF = MET_GetFieldRecord( "VParent", &fields );
    if( mF->defined )
      {
      tubeParentId = ( int )( mF->value[0] );
      }

    int tubeParentPoint = -1;
    mF = MET_GetFieldRecord( "Attachpt", &fields );
    if( mF->defined )
      {
      tubeParentPoint = ( int )( mF->value[0] );
      }

    float tubeColor[4];
    mF = MET_GetFieldRecord( "Color", &fields );
    tubeColor[0] = 1;
    tubeColor[1] = 0;
    tubeColor[2] = 0;
    tubeColor[3] = 1;
    if( mF->defined )
      {
      std::vector< float > colors;
      ::tube::StringToVector( ( char * )( mF->value ), colors, " " );
      for( unsigned int c=0; c<colors.size() && c<4; ++c )
        {
        tubeColor[ c ] = colors[c];
        }
      }

    char tubePointDim[80];
    strcpy( tubePointDim, "" );
    mF = MET_GetFieldRecord( "PointDim", &fields );
    if( mF->defined )
      {
      strcpy( tubePointDim, ( char * )( mF->value ) );
      }

    unsigned int nPoints = 0;
    mF = MET_GetFieldRecord( "NPoints", &fields );
    if( mF->defined )
      {
      nPoints = ( unsigned int )( mF->value[0] );
      }

    typename TubeType::Pointer tube = TubeType::New();

    bool useVoxelSize = true;
    for( unsigned int d=0; d<TDimension; ++d )
      {
      if( std::isnan(voxelSize[d]) || voxelSize[d] == 0 )
        {
        useVoxelSize = false;
        break;
        }
      }
    for( unsigned int j=0; j<nPoints; ++j )
      {
      typename TubeType::TubePointType pnt;

      typename TubeType::TubePointType::PointType x;
      for( unsigned int d=0; d<TDimension; ++d )
        {
        tmpReadStream >> x[d];
        tmpReadStream.get();
        if( useVoxelSize )
          {
          x[d] *= voxelSize[d];
          }
        }
      pnt.SetPositionInObjectSpace( x );

      double r;
      tmpReadStream >> r;
      pnt.SetRadiusInObjectSpace( r );

      tube->GetPoints().push_back( pnt );

      char c = ' ';
      while( c != '\n' && !tmpReadStream.eof() )
        {
        c = tmpReadStream.get();// to avoid unrecognize charactere
        }
      }

    tube->SetId( tubeId );
    tube->RemoveDuplicatePointsInObjectSpace();
    tube->ComputeTangentAndNormals();
    tube->GetProperty().SetColor( tubeColor[0], tubeColor[1],
      tubeColor[2] );
    tube->GetProperty().SetAlpha( tubeColor[3] );
    if( strlen( tubeAnat ) > 0 )
      {
      if( tubeAnat[0] == 'a' || tubeAnat[0] == 'A' )
        {
        tube->GetProperty().SetTagStringValue( "Artery", "True" );
        }
      else
        {
        tube->GetProperty().SetTagStringValue( "Artery", "False" );
        }
      }
    if( strlen( tubeTreeType ) > 0 )
      {
      if( tubeTreeType[0] == 'r' || tubeTreeType[0] == 'R' )
        {
        tube->SetRoot( true );
        }
      else
        {
        tube->SetRoot( false );
        if( tubeTreeType[0] == 'c' || tubeTreeType[0] == 'C' )
          {
          tube->SetParentId( tubeParentId );
          tube->SetParentPoint( tubeParentPoint );
          }
        }
      }

    m_TubeGroup->AddChild( tube );
    }

  tmpReadStream.close();

  for( unsigned int i=0; i<fields.size(); ++i )
    {
    delete fields[i];
    }
  fields.clear();

  return true;
}

template< unsigned int TDimension >
bool
TubeXIO< TDimension >
::Write( const std::string & _fileName )
{
  std::ofstream tmpWriteStream( _fileName.c_str(), std::ios::binary |
    std::ios::out );

  tmpWriteStream << "NDims: " << TDimension << std::endl;
  tmpWriteStream << "Dimensions: ";
  for( unsigned int i = 0; i < TDimension; ++i )
    {
    tmpWriteStream << this->m_Dimensions[i];
    if( i != TDimension - 1 )
      {
      tmpWriteStream << " ";
      }
    }
  tmpWriteStream << std::endl;

  char soType[80];
  snprintf( soType, 79, "Tube" );
  typename TubeType::ChildrenListType * tubeList =
    m_TubeGroup->GetChildren( 99999, soType );
  tmpWriteStream << std::endl;

  unsigned int nObjects = tubeList->size();
  tmpWriteStream << "NObjects: " << nObjects << std::endl;

  tmpWriteStream << "End Header:" << std::endl << std::endl;

  typename TubeType::ChildrenListType::iterator tIt = tubeList->begin();
  typename TubeType::ChildrenListType::iterator tEnd = tubeList->end();
  while( tIt != tEnd )
    {
    typename TubeType::Pointer tube = ( ( TubeType * )( tIt->GetPointer() ) );

    tmpWriteStream << "ID: " << tube->GetId() << std::endl;
    tmpWriteStream << "Type: Tube" << std::endl;
    std::string artery = "False";
    if( tube->GetProperty().GetTagStringValue( "Artery", artery ) &&
      artery == "True" )
      {
      tmpWriteStream << "Anat: artery" << std::endl;
      }
    if( tube->GetRoot() )
      {
      tmpWriteStream << "TreeType: root" << std::endl;
      }
    else if( tube->GetParentPoint() > 0 )
      {
      tmpWriteStream << "TreeType: child" << std::endl;
      tmpWriteStream << "VParent: " << tube->GetParentId() << std::endl;
      tmpWriteStream << "Attachpt: " << tube->GetParentPoint() << std::endl;
      }
    else
      {
      tmpWriteStream << "TreeType: orphan" << std::endl;
      }
    tmpWriteStream << "Color: "
      << tube->GetProperty().GetRed() << " "
      << tube->GetProperty().GetGreen() << " "
      << tube->GetProperty().GetBlue() << std::endl;
    tmpWriteStream << "PointDim: 4 x y z r" << std::endl;
    unsigned int nPoints = tube ->GetNumberOfPoints();
    tmpWriteStream << "NPoints: " << nPoints << std::endl;
    tmpWriteStream << "Points:" << std::endl;
    for( unsigned int i=0; i<nPoints; ++i )
      {
      typedef const typename TubeType::TubePointType  TubePointType;
      TubePointType * pnt;
      pnt = static_cast< TubePointType * >( tube->GetPoint( i ) );
      for( unsigned int k=0; k<TDimension; ++k )
        {
        tmpWriteStream << pnt->GetPositionInObjectSpace()[k] << " ";
        }
      tmpWriteStream << pnt->GetRadiusInObjectSpace() << std::endl;
      }
    tmpWriteStream << std::endl;

    ++tIt;
    }
  tmpWriteStream.close();

  tubeList->clear();
  delete tubeList;

  return true;
}

template< unsigned int TDimension >
void
TubeXIO< TDimension >
::SetTubeGroup( TubeGroupType * _tubes )
{
  m_TubeGroup = _tubes;
}

template< unsigned int TDimension >
typename GroupSpatialObject< TDimension >::Pointer &
TubeXIO< TDimension >
::GetTubeGroup( void )
{
  return m_TubeGroup;
}

} // tube namespace

} // itk namespace

#endif
