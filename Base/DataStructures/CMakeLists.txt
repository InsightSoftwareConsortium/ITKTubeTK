set(include_dirs
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR})

include_directories(${include_dirs})

if(BUILD_TESTING)
  add_subdirectory(Testing)
endif(BUILD_TESTING)
