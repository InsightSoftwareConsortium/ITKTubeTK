set( MODULE_NAME ClipTubes )
project( ${MODULE_NAME} )

if( NOT TubeTK_SOURCE_DIR )
  find_package( TubeTK REQUIRED )
  include( ${TubeTK_USE_FILE} )
endif( NOT TubeTK_SOURCE_DIR )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

find_package( ITK REQUIRED )
if( TubeTK_BUILD_WITHIN_SLICER )
  set( ITK_NO_IO_FACTORY_REGISTER_MANAGER 1 ) # See Libs/ITKFactoryRegistration/CMakeLists.txt
endif( TubeTK_BUILD_WITHIN_SLICER )
include( ${ITK_USE_FILE} )

SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES
    ${ITK_LIBRARIES} TubeTKIO TubeTKNumerics ITKIOMeta ITKIOSpatialObjects
  INCLUDE_DIRECTORIES
    ${TubeTK_SOURCE_DIR}/Base/CLI
    ${TubeTK_SOURCE_DIR}/Base/Common
    ${TubeTK_SOURCE_DIR}/Base/IO
    ${TubeTK_SOURCE_DIR}/Base/Filtering 
    ${TubeTK_SOURCE_DIR}/Base/Numerics )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
