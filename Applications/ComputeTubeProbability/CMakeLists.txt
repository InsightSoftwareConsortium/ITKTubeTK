project( tubeDensityProbability )
cmake_minimum_required( VERSION 2.8 )

# Find ITK
find_package( ITK REQUIRED )
include( ${USE_ITK_FILE} )

# Find GenerateCLP
find_package( GenerateCLP REQUIRED )
include( ${GenerateCLP_USE_FILE} )

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/ObjectDocuments
  )

set( PROJECT_SOURCE ${PROJECT_NAME}.cxx )

generateclp( PROJECT_SOURCE ${PROJECT_NAME}.xml )

add_executable( ${PROJECT_NAME} ${PROJECT_SOURCE} )
target_link_libraries( ${PROJECT_NAME} ${ITK_LIBRARIES} ObjectDocuments)
