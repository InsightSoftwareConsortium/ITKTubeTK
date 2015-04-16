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

#ifndef __itktubeTubeXIO_hxx
#define __itktubeTubeXIO_hxx

#include "itktubeTubeXIO.h"

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
  for ( unsigned int i = 0; i < TDimension; ++i )
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
  MET_InitReadField(mF, "NDims", MET_INT, false);
  fields.push_back(mF);
  int nDimsRecNum = MET_GetFieldRecordNumber("NDims", &fields);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "Dimensions", MET_FLOAT_ARRAY, false, nDimsRecNum);
  fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "VoxelSize", MET_FLOAT_ARRAY, false, nDimsRecNum);
  fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "NObjects", MET_INT, false);
  fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "End Header", MET_NONE, false);
  mF->terminateRead = true;
  fields.push_back(mF);

  std::ifstream tmpReadStream( _fileName.c_str(), std::ios::binary |
    std::ios::in );

  if( !tmpReadStream.rdbuf()->is_open() )
    {
    return false;
    }

  std::vector< MET_FieldRecordType * > extraFields;
  if( !MET_Read( tmpReadStream, &fields, ':', false, true, &extraFields ) )
    {
    std::cerr << "Tube: Read failed" << std::endl;
    return false;
    }

  mF = MET_GetFieldRecord( "NDims", &fields );
  if( mF->defined )
    {
    int nDims = (int)mF->value[0];

    if( (int)mF->value[0] != TDimension )
      {
      std::cerr << "Tube: Read failed: object is " << nDims
        << " dimensional and was expecting " << TDimension << " dimensional."
        << std::endl;
      return false;
      }
    }

  mF = MET_GetFieldRecord( "Dimensions", &fields );
  if(mF->defined)
    {
    for( unsigned int i=0; i<TDimension; ++i )
      {
      this->m_Dimensions[i] = (int)( mF->value[i] );
      }
    }

  double voxelSize[ TDimension ];
  mF = MET_GetFieldRecord( "VoxelSize", &fields );
  if(mF->defined)
    {
    for( unsigned int i=0; i<TDimension; ++i )
      {
      voxelSize[i] = (float)( mF->value[i] );
      }
    }

  int nObjects = 0;
  mF = MET_GetFieldRecord( "NObjects", &fields );
  if(mF->defined)
    {
    nObjects = (int)mF->value[0];
    }

  fields.clear();

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "ID", MET_INT, false);
  fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "Type", MET_STRING, false);
  fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "Anat", MET_STRING, false);
  fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "TreeType", MET_STRING, false);
  fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "Color", MET_STRING, false);
  fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "PointDim", MET_STRING, true);
  fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "NPoints", MET_INT, true);
  fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "Points", MET_NONE, true);
  mF->terminateRead = true;
  fields.push_back(mF);

  for( int obj = 0; obj < nObjects; ++obj )
    {
    if( !MET_Read( tmpReadStream, &fields, ':', false, true, &extraFields ) )
      {
      std::cerr << "Tube: Read failed" << std::endl;
      return false;
      }

    int tubeId = 0;
    mF = MET_GetFieldRecord( "ID", &fields );
    if(mF->defined)
      {
      tubeId = (int)( mF->value[0] );
      }

    char tubeType[80];
    strcpy( tubeType, "" );
    mF = MET_GetFieldRecord( "Type", &fields );
    if(mF->defined)
      {
      strcpy( tubeType, (char *)(mF->value) );
      }

    char tubeAnat[80];
    strcpy( tubeAnat, "" );
    mF = MET_GetFieldRecord( "Anat", &fields );
    if(mF->defined)
      {
      strcpy( tubeAnat, (char *)(mF->value) );
      }

    char tubeTreeType[80];
    strcpy( tubeTreeType, "" );
    mF = MET_GetFieldRecord( "TreeType", &fields );
    if(mF->defined)
      {
      strcpy( tubeTreeType, (char *)(mF->value) );
      }

    char tubeColor[80];
    strcpy( tubeColor, "" );
    mF = MET_GetFieldRecord( "Color", &fields );
    if(mF->defined)
      {
      strcpy( tubeColor, (char *)(mF->value) );
      }

    char tubePointDim[80];
    strcpy( tubePointDim, "" );
    mF = MET_GetFieldRecord( "PointDim", &fields );
    if(mF->defined)
      {
      strcpy( tubePointDim, (char *)(mF->value) );
      }

    unsigned int nPoints = 0;
    mF = MET_GetFieldRecord( "NPoints", &fields );
    if(mF->defined)
      {
      nPoints = (unsigned int)( mF->value[0] );
      }

    typename TubeType::Pointer tube = TubeType::New();

    for( unsigned int j=0; j<nPoints; ++j )
      {
      typename TubeType::TubePointType pnt;

      typename TubeType::TubePointType::PointType x;
      for( unsigned int d=0; d<TDimension; ++d )
        {
        tmpReadStream >> x[d];
        tmpReadStream.get();
        }
      pnt.SetPosition( x );

      double r;
      tmpReadStream >> r;
      pnt.SetRadius( r );

      tube->GetPoints().push_back( pnt );

      char c = ' ';
      while( c != '\n' && !tmpReadStream.eof() )
        {
        c = tmpReadStream.get();// to avoid unrecognize charactere
        }
      }

    tube->SetSpacing( voxelSize );
    tube->SetId( tubeId );
    tube->ComputeTangentAndNormals();

    m_TubeGroup->AddSpatialObject( tube );
    }

  tmpReadStream.close();

  return true;
}

template< unsigned int TDimension >
bool
TubeXIO< TDimension >
::Write( const std::string & _fileName )
{
  typedef typename TubeType::PointListType  PointListType;

  std::ofstream tmpWriteStream( _fileName.c_str(), std::ios::binary |
    std::ios::out);

  tmpWriteStream << "NDims: " << TDimension << std::endl;
  tmpWriteStream << "Dimensions: ";
  for ( unsigned int i = 0; i < TDimension; ++i )
    {
    tmpWriteStream << this->m_Dimensions[i];
    if ( i != TDimension - 1 )
      {
      tmpWriteStream << " ";
      }
    }
  tmpWriteStream << std::endl;

  char soType[80];
  sprintf( soType, "Tube" );
  typename TubeType::ChildrenListType * tubeList =
    m_TubeGroup->GetChildren( 99999, soType );
  tmpWriteStream << "VoxelSize:";
  for( unsigned int i=0; i<TDimension; ++i )
    {
    tmpWriteStream << " " << ( *(tubeList->begin()) )->GetSpacing()[i];
    }
  tmpWriteStream << std::endl;

  unsigned int nObjects = tubeList->size();
  tmpWriteStream << "NObjects: " << nObjects << std::endl;

  tmpWriteStream << "End Header:" << std::endl << std::endl;

  typename TubeType::ChildrenListType::iterator tIt = tubeList->begin();
  typename TubeType::ChildrenListType::iterator tEnd = tubeList->end();
  while( tIt != tEnd )
    {
    typename TubeType::Pointer tube = ((TubeType *)(tIt->GetPointer()));

    tmpWriteStream << "ID: " << tube->GetId() << std::endl;
    tmpWriteStream << "Type: Tube" << std::endl;
    tmpWriteStream << "Anat: artery" << std::endl;
    tmpWriteStream << "TreeType: orphan" << std::endl;
    tmpWriteStream << "Color: 1 0.3 0.21" << std::endl;
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
        tmpWriteStream << pnt->GetPosition()[k] << " ";
        }
      tmpWriteStream << pnt->GetRadius() << std::endl;
      }
    tmpWriteStream << std::endl;

    ++tIt;
    }
  tmpWriteStream.close();

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
