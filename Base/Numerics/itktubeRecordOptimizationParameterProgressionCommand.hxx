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

#ifndef __itktubeRecordOptimizationParameterProgressionCommand_hxx
#define __itktubeRecordOptimizationParameterProgressionCommand_hxx

#include "itktubeRecordOptimizationParameterProgressionCommand.h"

#include <itkSingleValuedNonLinearOptimizer.h>

namespace itk
{

namespace tube
{

template< unsigned int VNumberOfParameters, class ParametersValueType >
RecordOptimizationParameterProgressionCommand< VNumberOfParameters,
  ParametersValueType >
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
    HOFFSET( ParameterIterationType, Parameters ), parametersMemoryDataType );
  memoryDataType.insertMember( "CostFunctionValue",
    HOFFSET( ParameterIterationType, CostFunctionValue ), H5::PredType::NATIVE_DOUBLE );
  memoryDataType.insertMember( "CostFunctionDerivative",
    HOFFSET( ParameterIterationType, CostFunctionDerivative ), parametersMemoryDataType );
  this->H5ParameterIterationType = memoryDataType;
}


template< unsigned int VNumberOfParameters, class ParametersValueType >
void
RecordOptimizationParameterProgressionCommand< VNumberOfParameters,
  ParametersValueType >
::Execute( Object *caller, const EventObject & event )
{
  this->Execute( (const itk::Object *)caller, event);
}


template< unsigned int VNumberOfParameters, class ParametersValueType >
H5::CompType
RecordOptimizationParameterProgressionCommand< VNumberOfParameters,
  ParametersValueType >
::GetH5ParameterIterationType() const
{
  return this->H5ParameterIterationType;
}


template< unsigned int VNumberOfParameters, class ParametersValueType >
void
RecordOptimizationParameterProgressionCommand< VNumberOfParameters,
  ParametersValueType >
::Execute( const Object * object, const EventObject & event )
{
  if( StartEvent().CheckEvent( &event ) )
    {
    this->m_CurrentIteration = 0;

    const SingleValuedNonLinearOptimizer * optimizer =
      dynamic_cast< const SingleValuedNonLinearOptimizer * >( object );
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
      this->m_ParameterProgression[this->m_CurrentIteration].CostFunctionDerivative[ii] =
        derivative[ii];
      }

    return;
    }

  if( IterationEvent().CheckEvent( &event ) )
    {
    ++this->m_CurrentIteration;

    const SingleValuedNonLinearOptimizer * optimizer =
      dynamic_cast< const SingleValuedNonLinearOptimizer * >( object );
    if( object == NULL )
      {
      itkExceptionMacro(<< "Could not cast to itk::Optimizer.");
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
      this->m_ParameterProgression[this->m_CurrentIteration].CostFunctionDerivative[ii] =
        derivative[ii];
      }

    return;
    }

  if( EndEvent().CheckEvent( &event ) )
    {
    this->WriteParameterProgressionToFile();
    return;
    }
}


template< unsigned int VNumberOfParameters, class ParametersValueType >
void
RecordOptimizationParameterProgressionCommand< VNumberOfParameters,
  ParametersValueType >
::WriteParameterProgressionToFile( void ) const
{
  hsize_t dataspaceDimension[] = { this->m_ParameterProgression.size() };
  const int dataspaceRank = 1;

  H5::DataSpace dataSpace( dataspaceRank, dataspaceDimension );

  if( this->m_FileName.empty() )
    {
    itkExceptionMacro( << "FileName must be set." );
    }
  H5::H5File * file = new H5::H5File( this->m_FileName, H5F_ACC_TRUNC );

  H5::DataSet * dataset = new H5::DataSet( file->createDataSet( "OptimizationParameterProgression",
      this->H5ParameterIterationType, dataSpace ) );
  dataset->write( &(this->m_ParameterProgression[0]), this->H5ParameterIterationType );

  delete dataset;
  delete file;
}

} // End namespace itk

} // End namespace tube

#endif // End !defined(__itktubeRecordOptimizationParameterProgressionCommand_hxx)
