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
  execute_process( COMMAND "${PYTHON_EXECUTABLE}"
                   "-c"
                   "exec(\"import sys\\nif hasattr(sys,'real_prefix'):\\n  sys.exit(1)\")"
                   RESULT_VARIABLE TubeTKVirtualEnvSetup_InVirtualenv
                   ERROR_VARIABLE TubeTKVirtualEnvSetup_Error )
  if( TubeTKVirtualEnvSetup_InVirtualenv )
    message( WARNING "Configuring a virtual environment from a virtual environment\
 does not copy the packages from the original environment.")
  endif()
  # Check if virtualenv is available in the current environment (using
  # a virtualenv from a different environment could lead to errors
  # (e.g. python2 vs python3)
  execute_process( COMMAND "${PYTHON_EXECUTABLE}" -c
                   "import virtualenv;print virtualenv.__file__"
                   RESULT_VARIABLE TubeTKVirtualEnvSetup_Failed
                   OUTPUT_VARIABLE TubeTKVirtualEnvSetup_Output
                   ERROR_VARIABLE TubeTKVirtualEnvSetup_Error )
  if( TubeTKVirtualEnvSetup_Failed )
    message(FATAL_ERROR "virtualenv is required. Please install manually (e.g. \
pip install virtualenv")
  endif()
  get_filename_component( virtualenv_filename ${TubeTKVirtualEnvSetup_Output} NAME_WE)
  get_filename_component( virtualenv_extension ${TubeTKVirtualEnvSetup_Output} EXT)
  get_filename_component( virtualenv_directory ${TubeTKVirtualEnvSetup_Output} DIRECTORY)
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
    RESULT_VARIABLE TubeTKVirtualEnvSetup_Failed
    ERROR_VARIABLE TubeTKVirtualEnvSetup_Error )
  if( TubeTKVirtualEnvSetup_Failed )
    message( FATAL_ERROR ${TubeTKVirtualEnvSetup_Error} )
  endif()
else( NOT EXISTS "${activate_full}" )
  message( STATUS "Testing virtualenv found at ${PythonVirtualEnvDir}" )
endif( NOT EXISTS "${activate_full}" )

# Temporarily change the Python executable to the virtual environment one
set( PYTHON_TESTING_EXECUTABLE ${PythonVirtualEnvDir}/${PythonVirtualEnvScriptsDir}/python
    CACHE INTERNAL "Python to the python executable to use in tests." FORCE )

set( PYTHON_TESTING_MODULES )
if( ${TubeTK_USE_NUMPY} )
  set( PYTHON_TESTING_MODULES numpy )
endif( ${TubeTK_USE_NUMPY} )

if( TubeTK_USE_IPYTHON_NOTEBOOKS )
  list( APPEND PYTHON_TESTING_MODULES
    tornado
    pyzmq
    ipython
    jinja2
    matplotlib )
    # ipython[zmq]
endif( TubeTK_USE_IPYTHON_NOTEBOOKS )

if( TubeTK_USE_PYQTGRAPH )
  # Note: PyQt4 or PySide are required as is python-opengl, but these are not
  # pip installable.
  include(TubeTKCheckPythonLibrary)
  TubeTKCheckPythonLibrary(PyQt4)
  TubeTKCheckPythonLibrary(PySide)
  TubeTKCheckPythonLibrary(python-opengl)
  if( NOT PyQt4_FOUND OR NOT PySize_FOUND OR NOT python-opengl_FOUND )
    message(WARNING "With TubeTK_USE_PYQTGRAPH=ON and BUILD_TESTING=ON, PyQt4 or PySide\
 are required as is python-opengl, but these are not pip installable.")
  endif()

  list( APPEND PYTHON_TESTING_MODULES
    pyside
    pyqtgraph
    scipy
    tables )
endif( TubeTK_USE_PYQTGRAPH )

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
