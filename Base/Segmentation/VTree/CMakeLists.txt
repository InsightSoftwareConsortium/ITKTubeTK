PROJECT(TubeExtractorLib)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})
 
SET(TUBE_EXTRACTOR_SRCS
  itkOptBrent1D.cxx
  itkOptGoldenMean1D.cxx
  itkOptParabolicFit1D.cxx
  itkOptimizer1D.cxx
  itkOptimizerND.cxx
  itkSpline1D.cxx
  itkSplineApproximation1D.cxx
  itkSplineND.cxx
  vErf.cxx
  vRInt.cxx
  )

ADD_LIBRARY(TubeExtractor STATIC ${TUBE_EXTRACTOR_SRCS})
TARGET_LINK_LIBRARIES(TubeExtractor ${ITK_LIBRARIES})

#SUBDIRS(Testing)

