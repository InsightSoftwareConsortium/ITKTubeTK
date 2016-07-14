# Test that environment variables are set and correspond to this environment
option( TubeTK_TEST_ENVIRONMENT
  "Test that environment variables ITK_BUILD_DIR and TubeTK_BUILD_DIR are set and correspond to this project." ON )
if( TubeTK_TEST_ENVIRONMENT )
  add_test( NAME TestEnvironmentVariables
            COMMAND ${CMAKE_COMMAND} -DITK_BUILD_DIR:PATH=${ITK_DIR} -DTubeTK_BUILD_DIR:PATH=${PROJECT_BINARY_DIR}
                    -P ${PROJECT_SOURCE_DIR}/CMake/TubeTKTestEnvironmentVariables.cmake )
endif()
