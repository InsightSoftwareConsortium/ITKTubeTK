project(TubeTK)
#-----------------------------------------------------------------------------
# Require CMake 2.8 for the superbuild functionality
cmake_minimum_required(VERSION 2.8)

project(TubeTK)

#-----------------------------------------------------------------------------
# Superbuild is used by default
option(TubeTK_USE_SUPERBUILD
  "Build TubeTK and the projects it depends on via SuperBuild.cmake."
  ON)

if(TubeTK_USE_SUPERBUILD)
  option(TubeTK_SUPERBUILD_ITK "Build ITK via SuperBuild." ON)
  if(NOT TubeTK_SUPERBUILD_ITK)
    find_package(ITK)
    if(ITK_FOUND)
      include(${ITK_USE_FILE})
    else(ITK_FOUND)
      message(FATAL_ERROR 
        "Cannot build TubeTK without ITK.  Please set ITK_DIR.")
    endif(ITK_FOUND)
  endif(NOT TubeTK_SUPERBUILD_ITK)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
endif(TubeTK_USE_SUPERBUILD)

find_package(ITK)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
else(ITK_FOUND)
  message(FATAL_ERROR 
    "Cannot build TubeTK without ITK.  Please set ITK_DIR.")
endif(ITK_FOUND)

add_subdirectory(Base)
add_subdirectory(Applications)

install(FILES README.txt
  DESTINATION doc
)

