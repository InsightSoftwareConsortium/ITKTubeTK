include(ExternalProject)

set(base "${CMAKE_BINARY_DIR}")
set_property(DIRECTORY PROPERTY EP_BASE ${base})

set(shared ON) # use for BUILD_SHARED_LIBS on all subsequent projects
set(testing OFF) # use for BUILD_TESTING on all subsequent projects
set(build_type "Debug")
if(CMAKE_BUILD_TYPE)
  set(build_type "${CMAKE_BUILD_TYPE}")
endif() 

set( TubeTK_DEPENDS "" )

set(gen "${CMAKE_GENERATOR}")

##
## Check if sytem ITK or superbuild ITK
##
if(NOT USE_SYSTEM_ITK)

  ##
  ## Insight
  ##
  set(proj Insight)
  ExternalProject_Add(${proj}
    CVS_REPOSITORY 
      ":pserver:anonymous:insight@public.kitware.com:/cvsroot/Insight" 
    CVS_MODULE "Insight"
    SOURCE_DIR Insight
    BINARY_DIR Insight-Build
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
      -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
      -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      -DBUILD_SHARED_LIBS:BOOL=${shared}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DITK_USE_REVIEW:BOOL=ON
      -DITK_USE_OPTIMIZED_REGISTRATION_METHODS:BOOL=ON
    INSTALL_COMMAND ""
    )
  set( ITK_DIR "${base}/Insight-Build" )
  
endif(NOT USE_SYSTEM_ITK)


##
## TCLAP
##
set(proj tclap)
ExternalProject_Add(${proj}
  SVN_REPOSITORY 
    "http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel/tclap"
  SOURCE_DIR tclap
  BINARY_DIR tclap-Build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
    -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_TESTING:BOOL=OFF
  INSTALL_COMMAND ""
  )
set( TCLAP_DIR "${base}/tclap-Build" )


##
## ModuleDescriptionParser 
##
set(proj ModuleDescriptionParser)
if(NOT USE_SYSTEM_ITK)
  # Depends on ITK if ITK was build using superbuild
  set(ModuleDescriptionParser_DEPENDS "Insight")
else(NOT USE_SYSTEM_ITK)
  set(ModuleDescriptionParser_DEPENDS "")
endif(NOT USE_SYSTEM_ITK)
ExternalProject_Add(${proj}
  SVN_REPOSITORY 
    "http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel/ModuleDescriptionParser"
  SOURCE_DIR ModuleDescriptionParser
  BINARY_DIR ModuleDescriptionParser-Build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
    -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_TESTING:BOOL=OFF
    -DITK_DIR:PATH=${ITK_DIR}
  INSTALL_COMMAND ""
  DEPENDS ${ModuleDescriptionParser_DEPENDS}
  )
set( ModuleDescriptionParser_DIR "${base}/ModuleDescriptionParser-Build" )


##
## GenerateCLP 
##
set(proj GenerateCLP)
ExternalProject_Add(${proj}
  SVN_REPOSITORY 
    "http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel/GenerateCLP"
  SOURCE_DIR GenerateCLP
  BINARY_DIR GenerateCLP-Build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
    -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_TESTING:BOOL=OFF
    -DITK_DIR:PATH=${ITK_DIR}
    -DTCLAP_DIR:PATH=${TCLAP_DIR}
    -DModuleDescriptionParser_DIR:PATH=${ModuleDescriptionParser_DIR}
  INSTALL_COMMAND ""
  DEPENDS
    "tclap"
    "ModuleDescriptionParser"
  )
set( GenerateCLP_DIR "${base}/GenerateCLP-Build" )
set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "GenerateCLP" )


##
## OpenIGTLink 
##
if( TubeTK_USE_OpenIGTLink )
  set(proj OpenIGTLink)
  ExternalProject_Add(${proj}
    SVN_REPOSITORY "http://svn.na-mic.org/NAMICSandBox/trunk/OpenIGTLink"
    SOURCE_DIR OpenIGTLink
    BINARY_DIR OpenIGTLink-Build
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
      -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
      -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      -DBUILD_SHARED_LIBS:BOOL=${shared}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
    INSTALL_COMMAND ""
    )
  set( OpenIGTLink_DIR "${CMAKE_BINARY_DIR}/OpenIGTLink-Build" )
  set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "OpenIGTLink" )
else( TubeTK_USE_OpenIGTLink )
  set( OpenIGTLink_DIR "" )
endif( TubeTK_USE_OpenIGTLink )


##
## CTK 
##
if( TubeTK_USE_CTK )
  find_package( Git )
  set( QT_MIN_VERSION "4.6.0" )
  set( QT_OFFICIAL_VERSION "4.6" )
  set( QT_REQUIRED TRUE )
  find_package( Qt4 )
  if( NOT QT4_FOUND )
   MESSAGE(SEND_ERROR 
     "Qt ${QT_MIN_VERSION} or greater not found. Please check the QT_QMAKE_EXECUTABLE variable." )
  endif( NOT QT4_FOUND )
  set(proj CTK)
  ExternalProject_Add(${proj}
    GIT_REPOSITORY "http://github.com/commontk/CTK.git"
    SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj}
    BINARY_DIR ${proj}-Build
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      -DBUILD_SHARED_LIBS:BOOL=${shared}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DCTK_USE_GIT_PROTOCOL:BOOL=TRUE
      -DCTK_LIB_Widgets:BOOL=ON
      -DCTK_LIB_Visualization/VTK/Widgets:BOOL=OFF
      -DCTK_LIB_PluginFramework:BOOL=OFF
      -DCTK_PLUGIN_org.commontk.eventbus:BOOL=OFF
      -Dgit_EXECUTABLE:FILEPATH=${GIT_EXECUTABLE}
      -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
    INSTALL_COMMAND ""
    )
  set( CTK_DIR "${CMAKE_BINARY_DIR}/${proj}-Build" )
  set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "CTK" )
endif( TubeTK_USE_CTK )


##
## TubeTK - Normal Build
##
set(proj TubeTK)
ExternalProject_Add(${proj}
  DOWNLOAD_COMMAND ""
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR TubeTK-Build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DBUILDNAME:STRING=${BUILDNAME}
    -DSITE:STRING=${SITE}
    -DMAKECOMMAND:STRING=${MAKECOMMAND}
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
    -DBUILD_EXAMPLES:BOOL=${BUILD_EXAMPLES}
    -DBUILD_TESTING:BOOL=${BUILD_TESTING}
    -DBUILD_DOCUMENTATION:BOOL=${BUILD_DOCUMENTATION}
    -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
    -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DTubeTK_USE_SUPERBUILD:BOOL=FALSE
    -DTubeTK_USE_KWSTYLE:BOOL=${TubeTK_USE_KWSTYLE}
    -DITK_DIR:PATH=${ITK_DIR}
    -DGenerateCLP_DIR:PATH=${GenerateCLP_DIR}
    -DTubeTK_USE_OpenIGTLink:BOOL=${TubeTK_USE_OpenIGTLink}
    -DOpenIGTLink_DIR:PATH=${OpenIGTLink_DIR}
    -DTubeTK_USE_CTK:BOOL=${TubeTK_USE_CTK}
    -DCTK_DIR:PATH=${CTK_DIR}
  INSTALL_COMMAND ""
  DEPENDS
    ${TubeTK_DEPENDS}
)

