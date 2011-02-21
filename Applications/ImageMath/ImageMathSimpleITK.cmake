

#------------------------------------------------------------------------------
# Find required packages
#

# Find ITK
find_package(ITK REQUIRED)
if(ITK_FOUND)
  include(${ITK_USE_FILE})
endif()

# Find SimpleITK
find_package(SimpleITK REQUIRED)
if(SimpleITK_FOUND)
  include(${SimpleITK_USE_FILE})
endif()

# Cache the TubeTK directories
#set(TubeTK_SOURCE_DIR NOTFOUND CACHE PATH "Location of the TubeTK source")
#set(TubeTK_BUILD_DIR NOTFOUND CACHE PATH "Location of the TubeTK build")


#------------------------------------------------------------------------------
# Set up output directories
#

# Put all files in a bin directory to keep things tidy.
# Output directories.
if(NOT LIBRARY_OUTPUT_PATH)
  set(LIBRARY_OUTPUT_PATH ${ImageMath_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all libraries.")
endif()
if(NOT EXECUTABLE_OUTPUT_PATH)
  set(EXECUTABLE_OUTPUT_PATH ${ImageMath_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
endif()
mark_as_advanced(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)


#------------------------------------------------------------------------------
# Options for testing
#
option ( USE_TESTING "Build testing" ON )

if ( USE_TESTING )
  enable_testing()
  include(CTest)
  set(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
endif()


#------------------------------------------------------------------------------
# Go to subdirectories
#
subdirs(Code)
subdirs(ImageMath)

if (BUILD_TESTING)
  subdirs(Testing)
endif()

