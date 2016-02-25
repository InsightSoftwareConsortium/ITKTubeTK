cmake_minimum_required( VERSION 2.8.9 )
project( ITKTubeTK )
set(ITKTubeTK_LIBRARIES ITKTubeTK)

if( NOT ITK_SOURCE_DIR )

  find_package( ITK REQUIRED )
  list( APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR} )
  include( ITKModuleExternal )

else()

  set( TubeTK_INSTALL_DEVELOPMENT 1 )
  set( TubeTK_INSTALL_INCLUDE_DIR ${ITK_INSTALL_INCLUDE_DIR} )
  set( TubeTK_BUILD_TESTING 0 )

  # Include depends file to set ITKTubeTK_DEPENDS in current scope
  include( tubetk-depends.cmake )

  # Explicitly enable TubeTKBase libraries needed by ITKTubeTK
  foreach( name ${${PROJECT_NAME}_DEPENDS} )
    set( TubeTK_BUILD_BASE_${name} 1 )
  endforeach()

  # Build TubeTKBase
  add_subdirectory(
    ${CMAKE_CURRENT_SOURCE_DIR}/../../Base
    ${CMAKE_CURRENT_BINARY_DIR}/TubeTKBase-build
    )

  itk_module_impl()

endif()
