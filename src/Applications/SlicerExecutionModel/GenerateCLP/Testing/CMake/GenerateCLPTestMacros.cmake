
# --------------------------------------------------------------------------
# Helper macro
function(print_command_as_string command working_dir)
  if(PRINT_COMMAND)
    set(command_as_string)
    foreach(elem ${command})
      set(command_as_string "${command_as_string} \"${elem}\"")
    endforeach()
    message(STATUS "COMMAND:${command_as_string}")
    message(STATUS "WORKING_DIR: ${working_dir}")
  endif()
endfunction()

# --------------------------------------------------------------------------
# Macro used to generate CLIs to test GenerateCLP with
macro(GenerateCLP_TEST_PROJECT)
  include(CMakeParseArguments)

  set(options)
  set(oneValueArgs NAME)
  set(multiValueArgs)
  cmake_parse_arguments(TEST_CLI "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  #-----------------------------------------------------------------------------
  # Enable C++11
  #-----------------------------------------------------------------------------
  #####
  ##  Set the default target language standard properties
  if(NOT CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11) # Supported values are ``11``, ``14``, and ``17``.
  endif()
  if(NOT CMAKE_CXX_STANDARD_REQUIRED)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
  endif()
  if(NOT CMAKE_CXX_EXTENSIONS)
    set(CMAKE_CXX_EXTENSIONS OFF)
  endif()
  set(VALID_CXX_STANDARDS "11" "14" "17")
  if(NOT CMAKE_CXX_STANDARD IN_LIST VALID_CXX_STANDARDS )
     message(FATAL_ERROR "CMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD} not in know standards list\n ${VALID_CXX_STANDARDS}.")
  endif()

  project(${TEST_CLI_NAME} VERSION 0.0.1 LANGUAGES C CXX)

  #-----------------------------------------------------------------------------
  if(GenerateCLP_USE_JSONCPP)
    set(CMAKE_MODULE_PATH ${JsonCpp_CMAKE_MODULE_PATH} ${CMAKE_MODULE_PATH}) # Needed to locate FindJsonCpp.cmake
  endif()

  #-----------------------------------------------------------------------------
  find_package(GenerateCLP NO_MODULE REQUIRED)
  include(${GenerateCLP_USE_FILE})

  #-----------------------------------------------------------------------------
  if(GenerateCLP_USE_JSONCPP)
    find_package(JsonCpp REQUIRED)
    include_directories(${JsonCpp_INCLUDE_DIRS})
  endif()

  set(_additional_link_libraries)
  if(GenerateCLP_USE_JSONCPP)
    list(APPEND _additional_link_libraries ${JsonCpp_LIBRARIES})
  endif()
  if(GenerateCLP_USE_SERIALIZER)
    list(APPEND _additional_link_libraries ${ParameterSerializer_LIBRARIES})
  endif()

  #-----------------------------------------------------------------------------
  # Build
  #-----------------------------------------------------------------------------
  set(${PROJECT_NAME}_SOURCE ${PROJECT_NAME}.cxx)
  GENERATECLP(${PROJECT_NAME}_SOURCE ${PROJECT_NAME}.xml)
  sem_add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCE})
  if(_additional_link_libraries)
    target_link_libraries(${PROJECT_NAME} ${_additional_link_libraries})
  endif()

  #-----------------------------------------------------------------------------
  # Test
  #-----------------------------------------------------------------------------
  include(CTest)

  set(TEMP ${PROJECT_BINARY_DIR}/Testing/Temporary)
  file(MAKE_DIRECTORY ${TEMP})
endmacro()
