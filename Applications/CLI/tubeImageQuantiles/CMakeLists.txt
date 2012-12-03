project( tubeImageQuantiles )
cmake_minimum_required( VERSION 2.8 )

# Help finding BOOST
set(BOOST_ROOT  /software/boost_1_51_0/)
set(Boost_ADDITIONAL_VERSIONS "1.51")

find_package(Boost 1.51 REQUIRED)
if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
endif()

# Find ITK
find_package( ITK REQUIRED )
include( ${USE_ITK_FILE} )

# Find GenerateCLP
find_package( GenerateCLP REQUIRED )
include( ${GenerateCLP_USE_FILE} )


set(USER_CMAKE_CXX_FLAGS "-ftemplate-depth-100")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${USER_CMAKE_CXX_FLAGS}")

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  )

set( PROJECT_SOURCE ${PROJECT_NAME}.cxx )

generateclp( PROJECT_SOURCE ${PROJECT_NAME}.xml )

add_executable( ${PROJECT_NAME} ${PROJECT_SOURCE} )
target_link_libraries( ${PROJECT_NAME} ${ITK_LIBRARIES} )
