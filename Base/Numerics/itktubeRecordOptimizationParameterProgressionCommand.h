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

#ifndef __itktubeRecordOptimizationParameterProgressionCommand_h
#define __itktubeRecordOptimizationParameterProgressionCommand_h

#include <itkCommand.h>
#include <itkOptimizerParameters.h>

// Make sure to use this version to avoid accidentally getting VTK's internal
// HDF, for example.
#include <itk_H5Cpp.h>

namespace itk
{

namespace tube
{

/**
 * \class RecordOptimizationParameterProgressionCommand
 * Record progression of parameter values throughout an optimization process.
 *
 * This Command should observe an itk::SingleValuedNonLinearOptimizer.
 * */
template< unsigned int VNumberOfParameters, class TParametersValue = double >
class RecordOptimizationParameterProgressionCommand : public Command
{
public:
  typedef RecordOptimizationParameterProgressionCommand Self;
  typedef Command                                       Superclass;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;

  itkNewMacro( Self );

  itkStaticConstMacro( NumberOfParameters,
    unsigned int,
    VNumberOfParameters );
  typedef TParametersValue    ParametersValueType;
  typedef ParametersValueType CostFunctionValueType;
  typedef unsigned int        NumberOfIterationsType;

  typedef OptimizerParameters< ParametersValueType > FixedParametersType;

  /** Set/Get fixed parameters associated with the optimization. */
  virtual void SetFixedParameters( const FixedParametersType &
    fixedParameters );
  virtual const FixedParametersType & GetFixedParameters( void ) const;


  /** Set/Get the output file HDF5 file name. */
  itkSetStringMacro( FileName );
  itkGetStringMacro( FileName );

  /** Get the HDF5 DataType for a parameter iteration. */
  H5::CompType GetH5ParameterIterationType() const;

  /** Write the parameter progression to an HDF5 file.  This can be called
   * explicitly.  Alternatively, if this command observes the EndEvent on the
   * Optimizer, then it will be called at the end of optimization. */
  void WriteParameterProgressionToFile( void ) const;

  /** Convenience method to have this Command Observe the StartEvent,
   * IterationEvent, and EndEvent of the given Optimizer. */
  void Observe( Object * optimizer );

protected:
  RecordOptimizationParameterProgressionCommand( void );

  void Execute( Object *caller, const EventObject & event );

  void Execute( const Object * object, const EventObject & event );

  struct ParameterIteration
    {
    NumberOfIterationsType    Iteration;
    ParametersValueType       Parameters[NumberOfParameters];
    CostFunctionValueType     CostFunctionValue;
    CostFunctionValueType     CostFunctionDerivative[NumberOfParameters];
    };
  typedef ParameterIteration ParameterIterationType;

  typedef std::vector< ParameterIterationType > ParameterProgressionType;

private:
  // Purposely not implemented
  RecordOptimizationParameterProgressionCommand(
    const RecordOptimizationParameterProgressionCommand & );
  // Purposely not implemented
  void operator=( const RecordOptimizationParameterProgressionCommand & );

  NumberOfIterationsType   m_CurrentIteration;
  ParameterProgressionType m_ParameterProgression;
  FixedParametersType      m_FixedParameters;

  std::string              m_FileName;

  H5::CompType             m_H5ParameterIterationType;
}; // End class RecordOptimizationParameterProgressionCommand

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeRecordOptimizationParameterProgressionCommand.hxx"
#endif

// End !defined( __itktubeRecordOptimizationParameterProgressionCommand_h )
#endif
