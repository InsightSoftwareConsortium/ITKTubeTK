#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.13.4...3.14.2)
cmake_policy(VERSION 3.13.4...3.14.2)
include(${CMAKE_CURRENT_SOURCE_DIR}/../CMake/GenerateCLPTestMacros.cmake)

#-----------------------------------------------------------------------------
# Configure & Build
#-----------------------------------------------------------------------------
GenerateCLP_TEST_PROJECT(
  NAME CLPExample1
)

#-----------------------------------------------------------------------------
# Tests
#-----------------------------------------------------------------------------
if(GenerateCLP_USE_SERIALIZER)
  # Test normal command line exec + serialization at the end
  add_test(NAME GenerateCLPSerializer-${PROJECT_NAME}-Test1
    COMMAND $<TARGET_FILE:${PROJECT_NAME}>
    --spatialsamples 3324
    --selection 4
    --interpolation sinc
    --learningrate 0.002,0.001,0.0007,0.0002
    --translationscale 20
    --serialize ./GenerateCLPSerializer-${PROJECT_NAME}-Test1.json
    2
    Head.mha
    ProgrammingHead.mha
    ShrunkHead.mha
    WORKING_DIRECTORY ${TEMP}
    )

  # Test deserialization with no other args
  add_test(NAME GenerateCLPSerializer-${PROJECT_NAME}-Test2
    COMMAND $<TARGET_FILE:${PROJECT_NAME}>
    --deserialize ./GenerateCLPSerializer-${PROJECT_NAME}-Test1.json
    --serialize ./GenerateCLPSerializer-${PROJECT_NAME}-Test2.json
    WORKING_DIRECTORY ${TEMP}
    )
  set_tests_properties( GenerateCLPSerializer-${PROJECT_NAME}-Test2
    PROPERTIES DEPENDS GenerateCLPSerializer-${PROJECT_NAME}-Test1)

  # Test with command line (flagged) args overwriting JSON args
  add_test(NAME GenerateCLPSerializer-${PROJECT_NAME}-Test3
    COMMAND $<TARGET_FILE:${PROJECT_NAME}>
    --deserialize ./GenerateCLPSerializer-${PROJECT_NAME}-Test1.json
    --serialize ./GenerateCLPSerializer-${PROJECT_NAME}-Test3.json
    --spatialsamples 2112
    --learningrate 0.001,0.001,0.0005,0.0003
    2
    Head.mha
    ProgrammingHead.mha
    ShrunkHead.mha
    WORKING_DIRECTORY ${TEMP}
    )
  set_tests_properties( GenerateCLPSerializer-${PROJECT_NAME}-Test3
    PROPERTIES DEPENDS GenerateCLPSerializer-${PROJECT_NAME}-Test1)

  # Test with command line (positionnal) args overwriting JSON args
  add_test(NAME GenerateCLPSerializer-${PROJECT_NAME}-Test4
    COMMAND $<TARGET_FILE:${PROJECT_NAME}>
    --deserialize ./GenerateCLPSerializer-${PROJECT_NAME}-Test1.json
    --serialize ./GenerateCLPSerializer-${PROJECT_NAME}-Test4.json
    5
    OtherHead.mha
    OtherProgrammingHead.mha
    OtherShrunkHead.mha
    WORKING_DIRECTORY ${TEMP}
    )
  set_tests_properties( GenerateCLPSerializer-${PROJECT_NAME}-Test4
    PROPERTIES DEPENDS GenerateCLPSerializer-${PROJECT_NAME}-Test1)

  # Test argument flag collision
  add_test(NAME GenerateCLPSerializer-${PROJECT_NAME}-Test5
    COMMAND $<TARGET_FILE:${PROJECT_NAME}>
    --deserialize ./GenerateCLPSerializer-${PROJECT_NAME}-Test1.json
    --serialize ./GenerateCLPSerializer-${PROJECT_NAME}-Test5.json
    --translationscale 30
    WORKING_DIRECTORY ${TEMP}
    )
  set_tests_properties( GenerateCLPSerializer-${PROJECT_NAME}-Test5
    PROPERTIES DEPENDS GenerateCLPSerializer-${PROJECT_NAME}-Test1 WILL_FAIL ON)
endif()
