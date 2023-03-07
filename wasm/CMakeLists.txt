cmake_minimum_required(VERSION 3.16)
project(itk-tubetk-wasm)

set(CMAKE_CXX_STANDARD 17)

find_package(ITK REQUIRED
  COMPONENTS
    WebAssemblyInterface
    TubeTK
    )
include(${ITK_USE_FILE})

add_executable(crop-image crop-image.cxx)
target_link_libraries(crop-image PUBLIC ${ITK_LIBRARIES})
