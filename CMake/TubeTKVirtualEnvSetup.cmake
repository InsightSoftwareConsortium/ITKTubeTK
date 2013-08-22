find_package( PythonInterp REQUIRED )

# This sets the virtual environment directory
set( PythonVirtualenvHome ${TubeTK_BINARY_DIR}/Temporary/PythonVirtualenv )

get_filename_component( activate_full "${PythonVirtualenvHome}/bin/activate"
  ABSOLUTE )
if( NOT EXISTS "${activate_full}" )
  find_program( VIRTUALENV NAMES virtualenv virtualenv.py )
  if( VIRTUALENV )
    set( virtualenv_script "${VIRTUALENV}" )
  else( VIRTUALENV )
    set( virtualenv_script
      "@TubeTK_SOURCE_DIR@/Utilities/Python/virtualenv/virtualenv.py" )
  endif( VIRTUALENV )
  # Create the virtual environment.
  execute_process( COMMAND "${PYTHON_EXECUTABLE}"
    "${virtualenv_script}"
      "--python=${PYTHON_EXECUTABLE}"
      "--system-site-packages"
      "${PythonVirtualenvHome}"
    RESULT_VARIABLE failed
    ERROR_VARIABLE error )
  if( failed )
    message( ERROR ${error} )
  endif()
else( NOT EXISTS "${activate_full}" )
  message( STATUS "Testing virtualenv found at ${PythonVirtualenvHome}" )
endif( NOT EXISTS "${activate_full}" )

# Temporarily change the Python executable to the virtual environment one
if( WIN32 )
  set( PYTHON_TESTING_EXECUTABLE ${PythonVirtualenvHome}/Scripts/python
    CACHE INTERNAL "Python to the python executable to use in tests." FORCE )
else( )
  set( PYTHON_TESTING_EXECUTABLE ${PythonVirtualenvHome}/bin/python
    CACHE INTERNAL "Python to the python executable to use in tests." FORCE )
endif()

set( PYTHON_TESTING_MODULES
  numpy )

if( TubeTK_USE_NOTEBOOKS )
  list( APPEND PYTHON_TESTING_MODULES
    tornado
    pyzmq
    ipython[zmq]
    jinja2
    matplotlib )
endif( TubeTK_USE_NOTEBOOKS )

if( TubeTK_USE_PYQTGRAPH )
  # Note: PyQt4 or PySide are required as is python-opengl, but these are not
  # pip installable.
  list( APPEND PYTHON_TESTING_MODULES
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
