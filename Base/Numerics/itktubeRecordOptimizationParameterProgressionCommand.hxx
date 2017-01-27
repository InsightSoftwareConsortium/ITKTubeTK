/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#ifndef __itktubeRecordOptimizationParameterProgressionCommand_hxx
#define __itktubeRecordOptimizationParameterProgressionCommand_hxx

#include "itktubeRecordOptimizationParameterProgressionCommand.h"

#include <itkSingleValuedNonLinearOptimizer.h>

namespace itk
{

namespace tube
{

template< unsigned int VNumberOfParameters, class TParametersValue >
RecordOptimizationParameterProgressionCommand< VNumberOfParameters,
  TParametersValue >
::RecordOptimizationParameterProgressionCommand( void )
{
  H5::Exception::dontPrint();
  this->m_CurrentIteration = 0;

  H5::CompType memoryDataType( sizeof( ParameterIterationType ) );
  memoryDataType.insertMember( "Iteration",
    HOFFSET( ParameterIterationType, Iteration ),
    H5::PredType::NATIVE_UINT );
  hsize_t parametersDimension[] = { NumberOfParameters };
  H5::ArrayType parametersMemoryDataType( H5::PredType::NATIVE_DOUBLE,
    1, parametersDimension );
  memoryDataType.insertMember( "Parameters",
    HOFFSET( ParameterIterationType, Parameters ),
    parametersMemoryDataType );
  memoryDataType.insertMember( "CostFunctionValue",
    HOFFSET( ParameterIterationType, CostFunctionValue ),
    H5::PredType::NATIVE_DOUBLE );
  memoryDataType.insertMember( "CostFunctionDerivative",
    HOFFSET( ParameterIterationType, CostFunctionDerivative ),
    parametersMemoryDataType );
  this->m_H5ParameterIterationType = memoryDataType;
}


template< unsigned int VNumberOfParameters, class TParametersValue >
void
RecordOptimizationParameterProgressionCommand< VNumberOfParameters,
  TParametersValue >
::SetFixedParameters( const FixedParametersType & fixedParameters )
{
  if( this->m_FixedParameters != fixedParameters )
    {
    this->m_FixedParameters = fixedParameters;
    this->Modified();
    }
}


template< unsigned int VNumberOfParameters, class TParametersValue >
const typename RecordOptimizationParameterProgressionCommand<
  VNumberOfParameters, TParametersValue >::FixedParametersType &
RecordOptimizationParameterProgressionCommand< VNumberOfParameters,
  TParametersValue >
::GetFixedParameters() const
{
  return this->m_FixedParameters;
}


template< unsigned int VNumberOfParameters, class TParametersValue >
H5::CompType
RecordOptimizationParameterProgressionCommand< VNumberOfParameters,
  TParametersValue >
::GetH5ParameterIterationType() const
{
  return this->m_H5ParameterIterationType;
}


template< unsigned int VNumberOfParameters, class TParametersValue >
void
RecordOptimizationParameterProgressionCommand< VNumberOfParameters,
  TParametersValue >
::Execute( Object *caller, const EventObject & event )
{
  this->Execute( ( const itk::Object * )caller, event );
}


template< unsigned int VNumberOfParameters, class TParametersValue >
void
RecordOptimizationParameterProgressionCommand< VNumberOfParameters,
  TParametersValue >
::Execute( const Object * object, const EventObject & event )
{
  const SingleValuedNonLinearOptimizer * optimizer =
    dynamic_cast< const SingleValuedNonLinearOptimizer * >( object );
  if( StartEvent().CheckEvent( &event ) )
    {
    this->m_CurrentIteration = 0;

    if( object == NULL )
      {
      itkExceptionMacro( << "Could not cast to itk::Optimizer." );
      }

    this->m_ParameterProgression.resize( this->m_CurrentIteration + 1 );

    this->m_ParameterProgression[this->m_CurrentIteration].Iteration =
      this->m_CurrentIteration;

    const SingleValuedNonLinearOptimizer::ParametersType & parameters =
      optimizer->GetInitialPosition();
    for( unsigned int ii = 0; ii < NumberOfParameters; ++ii )
      {
      this->m_ParameterProgression[this->m_CurrentIteration].Parameters[ii] =
        parameters[ii];
      }

    this->m_ParameterProgression[this->m_CurrentIteration].CostFunctionValue =
      optimizer->GetValue( parameters );

    SingleValuedCostFunction::DerivativeType derivative;
    optimizer->GetCostFunction()->GetDerivative( parameters, derivative );
    for( unsigned int ii = 0; ii < NumberOfParameters; ++ii )
      {
      this->m_ParameterProgression[ this->m_CurrentIteration ].
        CostFunctionDerivative[ii] = derivative[ii];
      }

    return;
    }

  if( IterationEvent().CheckEvent( &event ) )
    {
    ++this->m_CurrentIteration;

    if( object == NULL )
      {
      itkExceptionMacro( << "Could not cast to itk::Optimizer." );
      }

    this->m_ParameterProgression.resize( this->m_CurrentIteration + 1 );

    this->m_ParameterProgression[this->m_CurrentIteration].Iteration =
      this->m_CurrentIteration;

    const SingleValuedNonLinearOptimizer::ParametersType & parameters =
      optimizer->GetCurrentPosition();
    for( unsigned int ii = 0; ii < NumberOfParameters; ++ii )
      {
      this->m_ParameterProgression[this->m_CurrentIteration].Parameters[ii] =
        parameters[ii];
      }

    this->m_ParameterProgression[this->m_CurrentIteration].CostFunctionValue =
      optimizer->GetValue( parameters );

    SingleValuedCostFunction::DerivativeType derivative;
    optimizer->GetCostFunction()->GetDerivative( parameters, derivative );
    for( unsigned int ii = 0; ii < NumberOfParameters; ++ii )
      {
      this->m_ParameterProgression[ this->m_CurrentIteration ].
        CostFunctionDerivative[ii] = derivative[ii];
      }

    return;
    }

  if( EndEvent().CheckEvent( &event ) )
    {
    this->WriteParameterProgressionToFile();
    return;
    }
}


template< unsigned int VNumberOfParameters, class TParametersValue >
void
RecordOptimizationParameterProgressionCommand< VNumberOfParameters,
  TParametersValue >
::WriteParameterProgressionToFile( void ) const
{
  if( this->m_FileName.empty() )
    {
    itkExceptionMacro( << "FileName must be set." );
    }
  H5::H5File * file = new H5::H5File( this->m_FileName, H5F_ACC_TRUNC );

  hsize_t iterationDataspaceDimension[] =
    {
    this->m_ParameterProgression.size()
    };
  const int dataspaceRank = 1;
  H5::DataSpace iterationDataSpace( dataspaceRank,
    iterationDataspaceDimension );
  H5::DataSet * iterationDataset = new H5::DataSet(
    file->createDataSet( "OptimizationParameterProgression",
    this->m_H5ParameterIterationType, iterationDataSpace ) );
  iterationDataset->write( &( this->m_ParameterProgression[0] ),
    this->m_H5ParameterIterationType );

  hsize_t fixedDataspaceDimension[] = { this->m_FixedParameters.Size() };
  H5::DataSpace fixedDataSpace( dataspaceRank, fixedDataspaceDimension );

  H5::DataSet * fixedDataset = new H5::DataSet( file->createDataSet(
    "FixedParameters", H5::PredType::NATIVE_DOUBLE, fixedDataSpace ) );
  if( fixedDataspaceDimension[0] > 0 )
    {
    fixedDataset->write( &( this->m_FixedParameters.GetElement( 0 ) ),
      H5::PredType::NATIVE_DOUBLE );
    }

  delete iterationDataset;
  delete fixedDataset;
  delete file;
}


template< unsigned int VNumberOfParameters, class TParametersValue >
void
RecordOptimizationParameterProgressionCommand< VNumberOfParameters,
  TParametersValue >
::Observe( Object * optimizer )
{
  optimizer->AddObserver( StartEvent(), this );
  optimizer->AddObserver( IterationEvent(), this );
  optimizer->AddObserver( EndEvent(), this );
}

} // End namespace tube

} // End namespace itk

// End !defined( __itktubeRecordOptimizationParameterProgressionCommand_hxx )
#endif
