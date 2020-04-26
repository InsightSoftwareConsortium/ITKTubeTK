itk_module_test()

set( TubeTK_DATA_ROOT ${TubeTK_SOURCE_DIR}/test/Baseline )

if( TubeTK_USE_ARRAYFIRE )
  find_package( ArrayFire REQUIRED )
endif()

if( TubeTK_USE_VTK )
  find_package( VTK REQUIRED )
endif()

add_subdirectory(CompareTools)
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
