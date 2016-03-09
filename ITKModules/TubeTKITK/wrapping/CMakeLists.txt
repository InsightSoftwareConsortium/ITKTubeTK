
# This macro will recursively get the <target> include directories
# and set the variable named <include_dir_varname> accordingly.
macro(_get_all_target_include_directories target include_dir_varname)

  get_target_property(_target_interface_link_libraries ${target} INTERFACE_LINK_LIBRARIES)
  foreach(_linked_target ${_target_interface_link_libraries})
    if( TARGET ${_linked_target} )
      get_target_property(_target_include_dirs ${_linked_target} INTERFACE_INCLUDE_DIRECTORIES)
      if(_target_include_dirs)
        set(_cleaned_target_include_dirs )
        foreach(_include_dir ${_target_include_dirs})
          string(REGEX MATCH "\\$<(BUILD|INSTALL)_INTERFACE:(.+)>" _is_genex ${_include_dir})
          if(NOT _is_genex STREQUAL "")
            set(_include_dir "${CMAKE_MATCH_2}")
          endif()
          list(APPEND _cleaned_target_include_dirs ${_include_dir})
        endforeach()
        list(APPEND ${include_dir_varname} ${_cleaned_target_include_dirs})
        _get_all_target_include_directories(${_linked_target} ${include_dir_varname})
      endif()
    endif()
  endforeach()
endmacro()


# Add target include directories as a directory property to satisfy CastXML.
#
# Note: This will NOT be needed when the ITK wrapping infrasctuture will
#       be updated to work with "Modern CMake" targets.
#
_get_all_target_include_directories( TubeTKITK _include_dirs )
include_directories( ${_include_dirs} )

itk_wrap_module( TubeTKITK )
itk_auto_load_submodules()
itk_end_wrap_module()
