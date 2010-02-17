project (HessianTubeness)

cmake_minimum_required(VERSION 2.8)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND CMAKE_POLICY)

# Disable MSVC 8 warnings
if(WIN32)
  option(DISABLE_MSVC8_DEPRECATED_WARNINGS "Disable Visual Studio 8 deprecated warnings" ON)
  mark_as_advanced(FORCE DISABLE_MSVC8_DEPRECATED_WARNINGS)  
  if(DISABLE_MSVC8_DEPRECATED_WARNINGS)
    add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
  endif(DISABLE_MSVC8_DEPRECATED_WARNINGS)
endif(WIN32)

# Find ITK
find_package(ITK REQUIRED)
include(${USE_ITK_FILE})

# Find GenerateCLP
find_package(GenerateCLP REQUIRED)
include(${GenerateCLP_USE_FILE})

include_directories(
  ${TubeTK_SOURCE_DIR}/Base/Preprocessing
  ${TubeTK_SOURCE_DIR}/Base
)

set ( HessianTubeness_SOURCE 
HessianTubeness.cxx
)

generateclp(HessianTubeness_SOURCE HessianTubeness.xml)

add_executable(HessianTubeness ${HessianTubeness_SOURCE})
target_link_libraries(HessianTubeness ${ITK_LIBRARIES})

install(TARGETS HessianTubeness
  RUNTIME DESTINATION bin
)
