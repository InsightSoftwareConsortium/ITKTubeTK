itk_module_test()

set( TubeTK_DATA_ROOT ${TubeTK_SOURCE_DIR}/test/Baseline )

add_subdirectory(CompareTools)
add_subdirectory(Common)
add_subdirectory(Filtering)
add_subdirectory(IO)
add_subdirectory(MetaIO)
add_subdirectory(Numerics)
#add_subdirectory(Python)
add_subdirectory(Registration)
add_subdirectory(Segmentation)
