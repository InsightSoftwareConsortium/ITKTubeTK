project( TubeTK )

#-----------------------------------------------------------------------------
# Require CMake 2.8 for the superbuild functionality
cmake_minimum_required( VERSION 2.8)

#-----------------------------------------------------------------------------
# Paths for TubeTK
INCLUDE_DIRECTORIES( 
  Base/Preprocessing
  )


#-----------------------------------------------------------------------------
# Dashboard setup
INCLUDE( ${CMAKE_ROOT}/Modules/Dart.cmake )
MARK_AS_ADVANCED( TCL_TCLSH DART_ROOT )
ENABLE_TESTING( )
IF( BUILD_TESTING )
  SET( BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard" )
  MARK_AS_ADVANCED( BUILDNAME )
ENDIF( BUILD_TESTING )
SET( BUILD_WARNING_REPORT_LIMIT -1)
SET( BUILD_ERROR_REPORT_LIMIT -1)
INCLUDE( "CTestConfig.cmake" )
CONFIGURE_FILE( ${TubeTK_SOURCE_DIR}/CTestCustom.cmake.in
                ${TubeTK_BINARY_DIR}/CTestCustom.cmake )


#-----------------------------------------------------------------------------
# Doxygen setup
add_subdirectory( Doxygen )

#-----------------------------------------------------------------------------
# KWStyle setup
OPTION( TubeTK_USE_KWSTYLE "Run KWStyle in order to check for violations of the coding standard." OFF )
MARK_AS_ADVANCED( TubeTK_USE_KWSTYLE )
IF( TubeTK_USE_KWSTYLE )
  FIND_PROGRAM( KWSTYLE_EXECUTABLE
    NAMES KWStyle 
    PATHS
    /usr/local/bin
    )
  MARK_AS_ADVANCED( KWSTYLE_EXECUTABLE )

  OPTION( KWSTYLE_USE_VIM_FORMAT "Set KWStyle to generate errors with a VIM-compatible format." OFF )
  OPTION( KWSTYLE_USE_MSVC_FORMAT "Set KWStyle to generate errors with a VisualStudio-compatible format." OFF )

  IF( KWSTYLE_USE_VIM_FORMAT )
    SET( KWSTYLE_ARGUMENTS -vim ${KWSTYLE_ARGUMENTS} )
  ENDIF( KWSTYLE_USE_VIM_FORMAT )

  IF( KWSTYLE_USE_MSVC_FORMAT )
    SET( KWSTYLE_ARGUMENTS -msvc ${KWSTYLE_ARGUMENTS} )
  ENDIF( KWSTYLE_USE_MSVC_FORMAT )

  CONFIGURE_FILE( ${TubeTK_SOURCE_DIR}/KWStyle/TubeTK.kws.xml.in
                  ${TubeTK_BINARY_DIR}/KWStyle-TubeTK.kws.xml )

  CONFIGURE_FILE( ${TubeTK_SOURCE_DIR}/KWStyle/TubeTKFiles.txt.in
                  ${TubeTK_BINARY_DIR}/KWStyle-TubeTKFiles.txt )

  CONFIGURE_FILE( ${TubeTK_SOURCE_DIR}/KWStyle/TubeTKOverwrite.txt.in
                  ${TubeTK_BINARY_DIR}/KWStyle-TubeTKOverwrite.txt )

  ADD_CUSTOM_TARGET( StyleCheck 
    COMMAND ${KWSTYLE_EXECUTABLE} -lesshtml -xml ${TubeTK_BINARY_DIR}/KWStyle-TubeTK.kws.xml -o ${TubeTK_BINARY_DIR}/KWStyle-TubeTKOverwrite.txt -dart ${TubeTK_BINARY_DIR} -1 1 -D ${TubeTK_BINARY_DIR}/KWStyle-TubeTKFiles.txt 
    WORKING_DIRECTORY ${TubeTK_BINARY_DIR}
    COMMENT "Coding Style Checker" VERBATIM )
ENDIF( TubeTK_USE_KWSTYLE )

#-----------------------------------------------------------------------------
# MSVC setup
IF( WIN32 )
  ADD_DEFINITIONS( -D_SCL_SECURE_NO_WARNINGS )
  IF( NOT CMAKE_FIRST_PASS_COMPLETED )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP" CACHE STRING "Flags used by all build types" FORCE )
  ENDIF( NOT CMAKE_FIRST_PASS_COMPLETED )
  OPTION( DISABLE_MSVC8_DEPRECATED_WARNINGS 
         "Disable Visual Studio 8 deprecated warnings" ON )
  MARK_AS_ADVANCED( FORCE DISABLE_MSVC8_DEPRECATED_WARNINGS )  
  IF( DISABLE_MSVC8_DEPRECATED_WARNINGS )
    ADD_DEFINITIONS( -D_CRT_SECURE_NO_DEPRECATE )
  ENDIF( DISABLE_MSVC8_DEPRECATED_WARNINGS )
ENDIF( WIN32 )

#-----------------------------------------------------------------------------
# Mac setup
IF( APPLE )
  IF( NOT CMAKE_FIRST_PASS_COMPLETED )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-common" CACHE STRING "Flags used by all build types" FORCE )
    SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-common" CACHE STRING "Flags used by all build types" FORCE )
  ENDIF( NOT CMAKE_FIRST_PASS_COMPLETED )
ENDIF( APPLE )

#-----------------------------------------------------------------------------
# Linux setup
IF( UNIX )
  IF( NOT CMAKE_FIRST_PASS_COMPLETED )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing" CACHE STRING "Flags used by all build types" FORCE ) 
    SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-strict-aliasing" CACHE STRING "Flags used by all build types" FORCE )
    IF( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
      SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC" CACHE STRING "Flags used by all build types" FORCE )
      SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC" CACHE STRING "Flags used by all build types" FORCE )
    ENDIF( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
  ENDIF( NOT CMAKE_FIRST_PASS_COMPLETED )
ENDIF( UNIX ) 

#-----------------------------------------------------------------------------
# Indicate that this cmakelists.txt file has already been parsed once.
# Keeps append commands from having an effect every time cmake is run.
SET( CMAKE_FIRST_PASS_COMPLETED ON CACHE BOOL "Modify the cache vars only on first pass thru cmake file" FORCE )
MARK_AS_ADVANCED( CMAKE_FIRST_PASS_COMPLETED )

#-----------------------------------------------------------------------------
# Superbuild setup
option( TubeTK_USE_SUPERBUILD
  "Build TubeTK and the projects it depends on via SuperBuild.cmake."
  ON )

if( TubeTK_USE_SUPERBUILD )

  option( USE_SYSTEM_ITK "Exclude ITK from SuperBuild and use an existing build instead." OFF )
  if( USE_SYSTEM_ITK )
    find_package( ITK REQUIRED )
    include(${ITK_USE_FILE})
  endif( USE_SYSTEM_ITK )
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()

else( TubeTK_USE_SUPERBUILD )

  find_package( ITK REQUIRED )
  include(${ITK_USE_FILE})

  find_package( GenerateCLP REQUIRED )
  include(${GenerateCLP_USE_FILE})

endif( TubeTK_USE_SUPERBUILD )

add_subdirectory( Base )
add_subdirectory( Applications )

install(
  FILES README.txt
  DESTINATION doc
)

