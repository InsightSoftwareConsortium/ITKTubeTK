include(ExternalData)

set( TubeTK_DATA_ROOT ${TubeTK_SOURCE_DIR}/data_keys)
set( TubeTKLib_DATA_ROOT ${TubeTKLib_SOURCE_DIR}/DataKeys)

if(NOT ExternalData_OBJECT_STORES)
  # Use ExternalData_OBJECT_STORES from environment as default.
  set(ExternalData_OBJECT_STORES_DEFAULT "")
  if(DEFINED "ENV{ExternalData_OBJECT_STORES}")
    file(TO_CMAKE_PATH "$ENV{ExternalData_OBJECT_STORES}"
      ExternalData_OBJECT_STORES_DEFAULT)
  endif()
endif()

# Select a data store.
if(NOT DEFINED ExternalData_OBJECT_STORES)
  if(DEFINED "ENV{ExternalData_OBJECT_STORES}")
    file(TO_CMAKE_PATH "$ENV{ExternalData_OBJECT_STORES}"
      ExternalData_OBJECT_STORES)
  else()
    if(DEFINED dashboard_data_name)
        set(ExternalData_OBJECT_STORES
          ${CTEST_DASHBOARD_ROOT}/${dashboard_data_name})
    else()
        set(ExternalData_OBJECT_STORES ${CTEST_DASHBOARD_ROOT}/ExternalData)
    endif()
  endif()
endif()

set(ExternalData_OBJECT_STORES "${ExternalData_OBJECT_STORES_DEFAULT}"
  CACHE STRING
  "Semicolon-separated list of data dirs in the layout %(algo)/%(hash).")
mark_as_advanced(ExternalData_OBJECT_STORES)
if(NOT ExternalData_OBJECT_STORES)
  set(ExternalData_OBJECT_STORES "${CMAKE_BINARY_DIR}/ExternalData/Objects")
  file(MAKE_DIRECTORY "${ExternalData_OBJECT_STORES}")
endif()
list(APPEND ExternalData_OBJECT_STORES
  "${CMAKE_SOURCE_DIR}/.ExternalData"
  )

set(ExternalData_BINARY_ROOT ${CMAKE_BINARY_DIR}/ExternalData)

# Expands %(algo:lower)
set(ExternalData_URL_TEMPLATES "" CACHE STRING
  "Additional URL templates for the ExternalData CMake script to look for testing data. E.g.
file:///var/bigharddrive/%(algo)/%(hash)")
mark_as_advanced(ExternalData_URL_TEMPLATES)
list(APPEND ExternalData_URL_TEMPLATES
  # Data published on Girder
  "https://data.kitware.com:443/api/v1/file/hashsum/%(algo)/%(hash)/download"
  )

# Tell ExternalData commands to transform raw files to content links.
# TODO: Condition this feature on presence of our pre-commit hook.
set(ExternalData_LINK_CONTENT SHA512)

# Emscripten currently has difficulty reading symlinks.
if(EMSCRIPTEN)
  set(ExternalData_NO_SYMLINKS 1)
endif()

# Match series of the form <base>.<ext>, <base>.<n>.<ext> such that <base> may
# end in a (test) number that is not part of any series numbering.
set(ExternalData_SERIES_PARSE "()(\\.[^./]*)$")
set(ExternalData_SERIES_MATCH "(\\.[0-9]+)?")

# Sometimes we want to download very large files.
set(ExternalData_TIMEOUT_ABSOLUTE 900)
