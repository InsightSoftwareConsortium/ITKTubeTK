set(proj ITK)

include(ExternalProject)

set(base "${CMAKE_BINARY_DIR}/CMakeExternals")
set_property(DIRECTORY PROPERTY EP_BASE ${base})

set(prefix "${base}/Install")

set(shared ON) # use for BUILD_SHARED_LIBS on all subsequent projects
set(testing OFF) # use for BUILD_TESTING on all subsequent projects
                                                                               set(build_type "")
if(CMAKE_BUILD_TYPE)
  set(build_type "${CMAKE_BUILD_TYPE}")
endif() 

ExternalProject_Add(${proj}
  LIST_SEPARATOR ${sep}
  CVS_REPOSITORY ":pserver:anonymous:insight@public.kitware.com:/cvsroot/Insight"                                                          
  CVS_MODULE "Insight"
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${prefix}
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_TESTING:BOOL=${testing}
    -DITK_USE_REVIEW:BOOL=ON
    -DITK_USE_OPTIMIZED_REGISTRATION_METHODS:BOOL=ON
)  

set(proj TubeTK-inner)

ExternalProject_Add(${proj}
  DOWNLOAD_COMMAND ""
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${prefix}
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -Ditk-based_USE_SUPERBUILD:BOOL=FALSE
    -DITK_DIR:PATH=${prefix}/lib/InsightToolkit
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR ${CMAKE_BINARY_DIR}/CMakeExternals/Build/TubeTK-inner
  DEPENDS
    "ITK"
)