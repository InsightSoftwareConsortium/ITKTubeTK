itk_module_test()

set( TubeTK_DATA_ROOT ${TubeTK_SOURCE_DIR}/test/Baseline )

add_subdirectory(CompareTools)
set( TubeTK_CompareTextFiles_EXE
  "$<TARGET_FILE:tubeCompareTextFilesTestDriver>" )
set( TubeTK_CompareImages_EXE
  "$<TARGET_FILE:tubeCompareImagesTestDriver>" )

add_subdirectory(Common)
add_subdirectory(Filtering)
add_subdirectory(IO)
add_subdirectory(MetaIO)
add_subdirectory(Numerics)
add_subdirectory(ObjectDocuments)
add_subdirectory(Registration)
add_subdirectory(Segmentation)

#if( ITK_WRAP_PYTHON )
  #add_subdirectory(Python)
#endif()
