project (HessianTubeness2D)

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
)

set ( HessianTubeness2D_SOURCE 
HessianTubeness2D.cxx
)

generateclp(HessianTubeness2D_SOURCE HessianTubeness2D.xml)

add_executable(HessianTubeness2D ${HessianTubeness2D_SOURCE})
target_link_libraries(HessianTubeness2D ${ITK_LIBRARIES})

install(TARGETS HessianTubeness2D
  RUNTIME DESTINATION bin
)
