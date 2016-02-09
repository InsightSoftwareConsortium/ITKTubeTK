cmake_minimum_required( VERSION 3.0.0 )

set(expected_nonempty_vars PROJECT_NAME GIT_REPOSITORY GIT_TAG MY_SOURCE_DIR)
foreach(varname IN LISTS expected_nonempty_vars)
  if("${${varname}}" STREQUAL "")
    message(FATAL_ERROR "Variable '${varname}' is empty")
  endif()
endforeach()

set(proj ${PROJECT_NAME})

project(DownloadSource LANGUAGES NONE)

include(ExternalProject)

ExternalProject_Add(${proj}
  GIT_REPOSITORY ${GIT_REPOSITORY}
  GIT_TAG ${GIT_TAG}
  SOURCE_DIR ${MY_SOURCE_DIR}
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  )
