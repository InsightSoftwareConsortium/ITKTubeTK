# Find LibSVM
#
# Sets LIBSVM_FOUND, LIBSVM_INCLUDE_DIR, LIBSVM_LIBRARY
#

if(LIBSVM_DIR)
  set(_svm_lib_dir "${LIBSVM_DIR}/lib")
  set(_svm_include_dir "${LIBSVM_DIR}/include")
endif()

find_path(LIBSVM_INCLUDE_DIR svm.h PATHS ${_svm_include_dir}  DOC "Path to the libsvm include directory")
find_library(LIBSVM_LIBRARY svm  PATHS ${_svm_lib_dir} DOC "The libsvm library.")

mark_as_advanced(LIBSVM_INCLUDE_DIR)
mark_as_advanced(LIBSVM_LIBRARY)

set(LIBSVM_FOUND 0)
if (LIBSVM_INCLUDE_DIR AND LIBSVM_LIBRARY)
  set(LIBSVM_FOUND 1)
  set(LIBSVM_LIBRARIES ${LIBSVM_LIBRARY})
endif()
