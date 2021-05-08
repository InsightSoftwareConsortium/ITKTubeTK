#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.13.4...3.14.2)
cmake_policy(VERSION 3.13.4...3.14.2)
include(${CMAKE_CURRENT_SOURCE_DIR}/../CMake/GenerateCLPTestMacros.cmake)

#-----------------------------------------------------------------------------
# Configure & Build
#-----------------------------------------------------------------------------
GenerateCLP_TEST_PROJECT(
  NAME CLPTestMultiple
)

#-----------------------------------------------------------------------------
# Tests
#-----------------------------------------------------------------------------
if(GenerateCLP_USE_SERIALIZER)
  list(APPEND FLAGS_TO_TEST
  --int
  --file
  --dir
  --image
  --geo
  --point
  --pointFile
  --roi
  )
  foreach(flag ${FLAGS_TO_TEST})
    add_test(NAME GenerateCLPSerializer-${PROJECT_NAME}${flag}-Test1
      COMMAND $<TARGET_FILE:${PROJECT_NAME}>
      ${flag} 1
      ${flag} 2,3,4
      --interruptionFlag
      ${flag} 5,6
      --serialize ./GenerateCLPSerializer-${PROJECT_NAME}${flag}-Test1.json
      WORKING_DIRECTORY ${TEMP}
      )

    # Test deserialization with no other args
    add_test(NAME GenerateCLPSerializer-${PROJECT_NAME}${flag}-Test2
      COMMAND $<TARGET_FILE:${PROJECT_NAME}>
      --deserialize ./GenerateCLPSerializer-${PROJECT_NAME}${flag}-Test1.json
      WORKING_DIRECTORY ${TEMP}
      )
    set_tests_properties( GenerateCLPSerializer-${PROJECT_NAME}${flag}-Test2
      PROPERTIES DEPENDS GenerateCLPSerializer-${PROJECT_NAME}${flag}-Test1)

    # Test deserialize with new arguments
    add_test(NAME GenerateCLPSerializer-${PROJECT_NAME}${flag}-Test3
      COMMAND $<TARGET_FILE:${PROJECT_NAME}>
      --testCommandLineOverWrite
      --deserialize ./GenerateCLPSerializer-${PROJECT_NAME}${flag}-Test1.json
      ${flag} -1,-2
      --interruptionFlag
      ${flag} -3
      ${flag} -4,-5,-6
      WORKING_DIRECTORY ${TEMP}
      )
    set_tests_properties( GenerateCLPSerializer-${PROJECT_NAME}${flag}-Test3
      PROPERTIES DEPENDS GenerateCLPSerializer-${PROJECT_NAME}${flag}-Test1)
  endforeach()

endif()
