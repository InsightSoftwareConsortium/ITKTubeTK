itk_module_test()

set( TubeTK_DATA_ROOT ${TubeTK_SOURCE_DIR}/test/Baseline )

#Waiting for SlicerExecutionModel to support ITKv5.1
#if( TubeTK_BUILD_APPLICATIONS )
  #add_subdirectory( Apps )
#endif()

add_subdirectory(CompareTools)
add_subdirectory(Common)
add_subdirectory(Filtering)
add_subdirectory(IO)
add_subdirectory(MetaIO)
add_subdirectory(Numerics)
add_subdirectory(ObjectDocuments)
#add_subdirectory(Python)
add_subdirectory(Registration)
add_subdirectory(Segmentation)
