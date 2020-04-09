project(TubeTKExamples)

find_package(ITK REQUIRED
  COMPONENTS TubeTK
    MinimalPathExtraction
    ITKIOImageBase
    ITKImageFunction
    ITKOptimizers
    ITKPath
    ITKIOMeta
    ITKIOPNG
    ITKIOJPEG
    )

include(${ITK_USE_FILE})

set( SlicerExecutionModel_DEFAULT_CLI_RUNTIME_OUTPUT_DIRECTORY
  "${ITK_BINARY_DIR}/bin" )
set( SlicerExecutionModel_DEFAULT_CLI_LIBRARY_OUTPUT_DIRECTORY
  "${ITK_BINARY_DIR}/lib" )
set( SlicerExecutionModel_DEFAULT_CLI_ARCHIVE_OUTPUT_DIRECTORY
  "${ITK_BINARY_DIR}/lib" )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )
mark_as_advanced( GENERATECLP_EXE )

add_subdirectory( Apps )
