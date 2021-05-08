# --------------------------------------------------------------------------
# Build tests
# --------------------------------------------------------------------------
set(KIT ${PROJECT_NAME})
create_test_sourcelist(Tests ${KIT}CxxTests.cxx
  ModuleDescriptionTest.cxx
  )

sem_add_executable(${KIT}CxxTests
  ModuleDescriptionTestingMacros.h
  ModuleDescriptionTestingUtilities.cxx
  ModuleDescriptionTestingUtilities.h
  ${Tests}
  )
target_link_libraries(${KIT}CxxTests ${KIT})

# --------------------------------------------------------------------------
# Build Parser1Test executable
# --------------------------------------------------------------------------
set(PARSER1TEST_SOURCE Parser1Test.cxx )
sem_add_executable(Parser1Test ${PARSER1TEST_SOURCE})
target_link_libraries(Parser1Test ModuleDescriptionParser )

# --------------------------------------------------------------------------
# Add tests
# --------------------------------------------------------------------------
set(TEST_DATA ${CMAKE_CURRENT_SOURCE_DIR}/TestData)
add_test(
  NAME Parser1Test1
  COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:Parser1Test> ${TEST_DATA}/ParserTest1.xml
  )

add_test(
  NAME Parser1Test2
  COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:Parser1Test> ${TEST_DATA}/ParserTest2.xml)

add_test(
  NAME Parser1Test3
  COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:Parser1Test> ${TEST_DATA}/ParserTest3.xml)

add_test(
  NAME ModuleDescriptionTest
  COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:${KIT}CxxTests> ModuleDescriptionTest ${TEST_DATA}
  )
