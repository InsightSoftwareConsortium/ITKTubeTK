
#
# Configure Test Prerequisites
#
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake/GenerateCLPTestPrerequisites.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMake/GenerateCLPTestPrerequisites.cmake
  @ONLY
  )

#
# Helper macro
#
set(_previous_test "NODEPENDS")
macro(generateclp_add_test cliname stepname)
  set(testname GenerateCLPTest-${cliname}-${stepname})
  add_test(NAME ${testname}
    COMMAND ${CMAKE_COMMAND}
    -DTEST_SOURCE_DIR:PATH=${CMAKE_CURRENT_SOURCE_DIR}/${cliname}
    -DTEST_BINARY_DIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/${cliname}
    -DTEST_CONFIGURATION:STRING=$<CONFIGURATION>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/CMake/GenerateCLPTest-${stepname}.cmake)
  if(NOT ${_previous_test} STREQUAL "NODEPENDS")
    set_tests_properties(${testname} PROPERTIES DEPENDS ${_previous_test})
  endif()
  set(_previous_test ${testname})
  set_property(TEST ${testname} PROPERTY LABELS ${PROJECT_NAME})
endmacro()

#
# Add Tests
#
set(_previous_test "NODEPENDS")
generateclp_add_test(Slicer3Example Configure)
generateclp_add_test(Slicer3Example Build)
generateclp_add_test(Slicer3Example Test)

set(_previous_test "NODEPENDS")
generateclp_add_test(PointExample1 Configure)
generateclp_add_test(PointExample1 Build)
#generateclp_add_test(PointExample1 Test)

set(_previous_test "NODEPENDS")
generateclp_add_test(CLPExample1 Configure)
generateclp_add_test(CLPExample1 Build)
generateclp_add_test(CLPExample1 Test)

set(_previous_test "NODEPENDS")
generateclp_add_test(CLPTestMultiple Configure)
generateclp_add_test(CLPTestMultiple Build)
generateclp_add_test(CLPTestMultiple Test)

if(SlicerExecutionModel_USE_UTF8)
  set(_previous_test "NODEPENDS")
  generateclp_add_test(Utf8Example Configure)
  generateclp_add_test(Utf8Example Build)
  generateclp_add_test(Utf8Example Test)
endif()
