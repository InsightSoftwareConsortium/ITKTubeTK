PROJECT(QtSlicer)

INCLUDE (${CMAKE_ROOT}/Modules/FindOpenGL.cmake)

#
#  Find where ITK was Built or Installed.
#
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)


FIND_PACKAGE(Qt4 REQUIRED) 

SET(QT_USE_QTOPENGL TRUE)
INCLUDE(${QT_USE_FILE})

MESSAGE(${QT_INCLUDE_DIR})

INCLUDE_DIRECTORIES(
 ${QT_INCLUDE_DIR}
 ${QT_INCLUDE_DIR}/Qt
 ${QtSlicer_BINARY_DIR}
 ${QtSlicer_SOURCE_DIR}
)


LINK_LIBRARIES(
 ${ITK_LIBRARIES}
 ${QT_QT_LIBRARY}
 ${QT_GL_LIBRARY}
 ${OPENGL_glu_LIBRARY}
 ${OPENGL_LIBRARY}
)

SET(QtSlicer_SRCS
  QtSlicerTest.cxx
  QtGlSliceView.cxx
  QtSlicer.cxx
)

QT4_WRAP_UI(QtSlicer_UIS_H
  QtSlicerGUI.ui
  QtSlicerHelpGUI.ui
  )

QT4_WRAP_CPP(QtSlicer_MOC_SRCS
  QtGlSliceView.h
  )

source_group("Generated" FILES
  ${QtSlicer_UIS_H}
  ${QtSlicer_MOC_SRCS}
  )

ADD_DEFINITIONS(-DQT_DLL)

ADD_DEFINITIONS(-DQT_THREAD_SUPPORT)

ADD_EXECUTABLE(${PROJECT_NAME} 
  ${QtSlicer_SRCS}
  ${QtSlicer_UIS_H}
  ${QtSlicer_MOC_SRCS}
  )

TARGET_LINK_LIBRARIES(${PROJECT_NAME}
  ${QT_LIBRARIES}
  )

slicer3_set_plugins_output_path( ${PROJECT_NAME} )
slicer3_install_plugins( ${PROJECT_NAME} )