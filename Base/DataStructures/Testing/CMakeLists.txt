set(KIT DataStructures)

create_test_sourcelist(Tests ${KIT}CxxTests.cxx
  itkGeometryObjectTest1.cxx)

add_executable(${KIT}CxxTests ${Tests})
target_link_libraries(${KIT}CxxTests ITKCommon)

add_test(itkGeometryObjectTest1 
  ${KIT}CxxTests itkGeometryObjectTest1)