# Create the fooCLP.h file.  This file should be included in
# foo.cxx and the PARSE_ARGS pre-processor macro used to process command line
# arguments.
#
# Usage: GENERATECLP(foo_SRCS XML_FILE [LOGO_FILE])
#
# foo_SRCS should be a list of source containing foo.cxx.  The XML_FILE should
# should have the same basename as the .cxx file, e.g. foo.xml.
macro(GENERATECLP SOURCES XML_FILE)

  # Sanity checks
  if("${GENERATECLP_EXE}" STREQUAL "")
    message(FATAL_ERROR "GENERATECLP_EXE value is invalid. It shouldn't be an empty string !")
  endif()

  # The shell into which nmake.exe executes the custom command has some issues
  # with mixing quoted and unquoted arguments :( Let's help.

  if(CMAKE_GENERATOR MATCHES "NMake Makefiles")
    set(verbatim "")
    set(quote "\"")
  else()
    set(verbatim "VERBATIM")
    set(quote "")
  endif()

  # What is the filename without the extension
  get_filename_component(TMP_FILENAME ${ARGV1} NAME_WE)

  # The input file might be full path so handle that
  get_filename_component(TMP_FILEPATH ${ARGV1} PATH)

  # Compute the input filename
  if(TMP_FILEPATH)
    set(TMP_INPUT ${TMP_FILEPATH}/${TMP_FILENAME}.xml)
  else()
    set(TMP_INPUT ${CMAKE_CURRENT_SOURCE_DIR}/${TMP_FILENAME}.xml)
  endif()

  # add custom command to output
  if ("x${ARGV2}" STREQUAL "x")

    ## CASE OF NO LOGO FILE
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}CLP.h
      DEPENDS "${GENERATECLP_EXE}" ${TMP_INPUT}
      COMMAND "${GENERATECLP_EXE}"
      --InputXML "${quote}${TMP_INPUT}${quote}"
      --OutputCxx "${quote}${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}CLP.h${quote}"
      ${verbatim}
      )

    foreach(INPUT_SOURCE_FILE ${${SOURCES}} ${TMP_FILENAME}.cxx)
      set_source_files_properties(${INPUT_SOURCE_FILE} PROPERTIES
        OBJECT_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}CLP.h)
    endforeach()

  else()

    ## CASE WITH LOGO FILE
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}CLP.h
      DEPENDS "${GENERATECLP_EXE}" ${TMP_INPUT} ${ARGV2}
      COMMAND "${GENERATECLP_EXE}"
      --logoFiles "${quote}${ARGV2}${quote}"
      --InputXML "${quote}${TMP_INPUT}${quote}"
      --OutputCxx "${quote}${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}CLP.h${quote}"
      ${verbatim}
      )

    foreach(INPUT_SOURCE_FILE ${${SOURCES}} ${TMP_FILENAME}.cxx)
      set_source_files_properties(
        ${INPUT_SOURCE_FILE} PROPERTIES
          OBJECT_DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}CLP.h
          OBJECT_DEPENDS ${ARGV2})
    endforeach()

    # Mark the logo include file as a header file
    set_source_files_properties(
      ${ARGV2}
      PROPERTIES HEADER_FILE_ONLY TRUE)

  endif()

  set_source_files_properties(
    ${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}CLP.h PROPERTIES
      HEADER_FILE_ONLY TRUE # Mark the file as a header file
      GENERATED TRUE # Mark the file as generated
      )

  set(${SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/${TMP_FILENAME}CLP.h ${${SOURCES}})
  include_directories(${CMAKE_CURRENT_BINARY_DIR})

endmacro()
