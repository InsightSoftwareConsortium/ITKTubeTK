cmake_minimum_required(VERSION 3.8)

set(P_NAME GenerateInputMipPointImageTest)

project(${P_NAME})

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

set(CMAKE_CXX_STANDARD 11)

add_executable(${P_NAME} ${P_NAME}.cxx)

target_link_libraries(${P_NAME} ${ITK_LIBRARIES})
