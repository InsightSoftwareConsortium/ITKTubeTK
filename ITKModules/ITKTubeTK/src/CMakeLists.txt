
set( ITKTubeTK_SRC
  itkTubeTK_fake_library_src.cxx
)

add_library( ITKTubeTK STATIC ${ITKTubeTK_SRC} )
target_link_libraries( ITKTubeTK PUBLIC TubeTKFiltering )

itk_module_link_dependencies()
itk_module_target( ITKTubeTK )
