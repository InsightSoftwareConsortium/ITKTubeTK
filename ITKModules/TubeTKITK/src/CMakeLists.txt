
set( TubeTKITK_SRC
  itkTubeTK_fake_library_src.cxx
)

# Include depends file to set TubeTKITK_LIBRARIES in current scope
include( ../tubetk-depends.cmake )

add_library( TubeTKITK STATIC ${TubeTKITK_SRC} )
target_link_libraries( TubeTKITK PUBLIC ${TubeTKITK_LIBRARIES} )

itk_module_link_dependencies()
itk_module_target( TubeTKITK )
