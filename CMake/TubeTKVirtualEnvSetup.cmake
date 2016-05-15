find_package( PythonInterp REQUIRED )

find_package( PythonLibs REQUIRED )

# This sets the virtual environment directory
if( NOT PythonVirtualEnvDir )
  set( PythonVirtualEnvDir ${TubeTK_BINARY_DIR}/Temporary/PythonVirtualenv )
endif()

if( WIN32 )
  set( PythonVirtualEnvScriptsDir Scripts)
else( WIN32 )
  set( PythonVirtualEnvScriptsDir bin)
endif( WIN32 )
get_filename_component( activate_full
  "${PythonVirtualEnvDir}/${PythonVirtualEnvScriptsDir}/activate" ABSOLUTE )

if( NOT EXISTS "${activate_full}" )

  # Check if we are in a virtal environment
  execute_process( COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "import sys; sys.exit(1) if hasattr(sys,'real_prefix') else sys.exit(0)"
    RESULT_VARIABLE inVirtualEnv
    ERROR_VARIABLE inVirtualEnv_error )
  if( inVirtualEnv )
    message( WARNING "Configuring a virtual environment from a virtual\
      environment does not copy the packages from the original environment.")
  endif()
  # Check if virtualenv is available in the current environment (using
  # a virtualenv from a different environment could lead to errors
  # (e.g. python2 vs python3)
  execute_process( COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "import virtualenv;print(virtualenv.__file__)"
    RESULT_VARIABLE virtualEnv_failed
    OUTPUT_VARIABLE virtualEnv_output
    ERROR_VARIABLE virtualEnv_error )
  if( virtualEnv_failed )
    message(FATAL_ERROR "virtualenv is required, e.g., pip install virtualenv")
  endif()
  get_filename_component( virtualenv_filename ${virtualEnv_output} NAME_WE)
  get_filename_component( virtualenv_extension ${virtualEnv_output} EXT)
  get_filename_component( virtualenv_directory ${virtualEnv_output} DIRECTORY)
  set( virtualenv_script ${virtualenv_directory}/${virtualenv_filename})
  if( virtualenv_extension )
    set( virtualenv_script ${virtualenv_script}.py )
    if (NOT EXISTS ${virtualenv_script} )
      message(FATAL_ERROR "${virtualenv_script} not found")
    endif()
  endif()
  message(STATUS "Virtualenv path - ${virtualenv_script}")
  # Create the virtual environment.
  execute_process( COMMAND "${PYTHON_EXECUTABLE}"
    "${virtualenv_script}"
      "--python=${PYTHON_EXECUTABLE}"
      "--system-site-packages"
      "${PythonVirtualEnvDir}"
    RESULT_VARIABLE createVirtualEnv_failed
    ERROR_VARIABLE createVirtualEnv_error )
  if( createVirtualEnv_failed )
    message( FATAL_ERROR ${createVirtualEnv_error} )
  endif()
else( NOT EXISTS "${activate_full}" )
  message( STATUS "Testing virtualenv found at ${PythonVirtualEnvDir}" )
endif( NOT EXISTS "${activate_full}" )

# Temporarily change the Python executable to the virtual environment one
set( PYTHON_TESTING_EXECUTABLE
  ${PythonVirtualEnvDir}/${PythonVirtualEnvScriptsDir}/python
  CACHE INTERNAL "Python to the python executable to use in tests." FORCE )

set( PYTHON_TESTING_MODULES )

# numpy
if( ${TubeTK_USE_NUMPY} )

  list( APPEND PYTHON_TESTING_MODULES
    scipy
    numpy )

endif( ${TubeTK_USE_NUMPY} )

# pyqtgraph
if( TubeTK_USE_PYQTGRAPH )

  list( APPEND PYTHON_TESTING_MODULES
    pyqtgraph )

endif( TubeTK_USE_PYQTGRAPH )

# ipython and other things used by examples
if( ${TubeTK_USE_PYTHON_EXAMPLES_AS_TESTS} )

  list( APPEND PYTHON_TESTING_MODULES
    ipython
    tornado
    pyzmq
    jinja2
    tables
    matplotlib )

endif( ${TubeTK_USE_PYTHON_EXAMPLES_AS_TESTS} )

configure_file(
  "${TubeTK_SOURCE_DIR}/CMake/PythonVirtualEnvInstall.cmake.in"
  "${TubeTK_BINARY_DIR}/CMake/PythonVirtualEnvInstall.cmake"
  @ONLY )

# Install the virtual environment
add_custom_target( PythonVirtualenvInstall ALL
  COMMAND ${CMAKE_COMMAND} -P
  "${TubeTK_BINARY_DIR}/CMake/PythonVirtualEnvInstall.cmake" )

# Setup IPython notebook driver
set( NOTEBOOK_TEST_DRIVER
  ${TubeTK_SOURCE_DIR}/Utilities/Python/EvaluateIPythonNotebook.py
  CACHE INTERNAL "Test driver command for IPython Notebooks." FORCE )
