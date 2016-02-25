
set( ITKTubeTK_SRC
  itkTubeTK_fake_library_src.cxx
)

# Include depends file to set ITKTubeTK_LIBRARIES in current scope
include( ../tubetk-depends.cmake )

add_library( ITKTubeTK STATIC ${ITKTubeTK_SRC} )
target_link_libraries( ITKTubeTK PUBLIC ${ITKTubeTK_LIBRARIES} )

itk_module_link_dependencies()
itk_module_target( ITKTubeTK )
