project( tubeLabelTransfer )
cmake_minimum_required( VERSION 2.8 )

# Find ITK
find_package( ITK REQUIRED )
include( ${USE_ITK_FILE} )

# Find GenerateCLP
find_package( GenerateCLP REQUIRED )
include( ${GenerateCLP_USE_FILE} )

# Release flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -msse3" ${CMAKE_FLAGS_RELEASE})

# Find Boost
#set(Boost_ADDITIONAL_VERSIONS "1.51")
#find_package(Boost 1.51 REQUIRED)
#include_directories(${Boost_INCLUDE_DIRS})

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/ObjectDocuments
  )

set( PROJECT_SOURCE ${PROJECT_NAME}.cxx)

generateclp( PROJECT_SOURCE ${PROJECT_NAME}.xml )

add_executable( ${PROJECT_NAME} ${PROJECT_SOURCE} )
target_link_libraries( ${PROJECT_NAME} ${ITK_LIBRARIES} ObjectDocuments)
