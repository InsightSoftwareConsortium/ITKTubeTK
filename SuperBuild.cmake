include(ExternalProject)

set(base "${CMAKE_BINARY_DIR}/CMakeExternals")
set_property(DIRECTORY PROPERTY EP_BASE ${base})

set(prefix "${base}/Install")

set(shared ON) # use for BUILD_SHARED_LIBS on all subsequent projects
set(testing OFF) # use for BUILD_TESTING on all subsequent projects
set(build_type "Debug")
if(CMAKE_BUILD_TYPE)
  set(build_type "${CMAKE_BUILD_TYPE}")
endif() 

if(NOT USE_SYSTEM_ITK)

set(proj ITK)

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
    -DBUILD_TESTING:BOOL=OFF
    -DITK_USE_REVIEW:BOOL=ON
    -DITK_USE_OPTIMIZED_REGISTRATION_METHODS:BOOL=ON
)

SET( ITK_DIR ${prefix}/lib/InsightToolkit )

endif(NOT USE_SYSTEM_ITK)

set(proj tclap)

ExternalProject_Add(${proj}
  SVN_REPOSITORY "http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel/tclap"
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${prefix}/tclapsandbox
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_TESTING:BOOL=OFF
)

set(proj ModuleDescriptionParser)
 
if(NOT USE_SYSTEM_ITK)
  ExternalProject_Add(${proj}
    SVN_REPOSITORY 
      "http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel/ModuleDescriptionParser"
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX:PATH=${prefix}
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_SHARED_LIBS:BOOL=${shared}
      -DBUILD_TESTING:BOOL=OFF
      -DITK_DIR:PATH=${ITK_DIR}
    DEPENDS
      "ITK"
  )
else(NOT USE_SYSTEM_ITK)
  ExternalProject_Add(${proj}
    SVN_REPOSITORY 
      "http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel/ModuleDescriptionParser"
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX:PATH=${prefix}
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_SHARED_LIBS:BOOL=${shared}
      -DBUILD_TESTING:BOOL=OFF
      -DITK_DIR:PATH=${ITK_DIR}
  )
endif(NOT USE_SYSTEM_ITK)

set(proj GenerateCLP)
 
ExternalProject_Add(${proj}
  SVN_REPOSITORY 
    "http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel/GenerateCLP"
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${prefix}
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_TESTING:BOOL=OFF
    -DITK_DIR:PATH=${ITK_DIR}
    -DTCLAP_DIR:PATH=${prefix}/tclapsandbox/lib/tclap
    -DModuleDescriptionParser_DIR:PATH=${prefix}/lib/ModuleDescriptionParser
  DEPENDS
    "tclap"
    "ModuleDescriptionParser"
)

set(proj TubeTK-inner)

ExternalProject_Add(${proj}
  DOWNLOAD_COMMAND ""
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${prefix}
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DTubeTK_USE_SUPERBUILD:BOOL=FALSE
    -DITK_DIR:PATH=${ITK_DIR}
    -DGenerateCLP_DIR:PATH=${prefix}/lib/GenerateCLP
    -DBUILD_TESTING:BOOL=${BUILD_TESTING}
    -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
    -DBUILD_DOCUMENTATION:BOOL=${BUILD_DOCUMENTATION}
    -DBUILDNAME:STRING=${BUILDNAME}
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
    -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
    -DMAKECOMMAND:STRING=${MAKECOMMAND}
    -DSITE:STRING=${SITE}
    -DTubeTK_USE_KWSTYLE:BOOL=${TubeTK_USE_KWSTYLE}
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR ${CMAKE_BINARY_DIR}/CMakeExternals/Build/TubeTK-inner
  DEPENDS
    "GenerateCLP"
)
