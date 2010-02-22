set(TubeTK_Base_TubeSegmentation_H_Files
  tubeCurves2DImageFilter.h
  )

set(TubeTK_Base_TubeSegmentation_TXX_Files
  )

set(TubeTK_Base_TubeSegmentation_CXX_Files
  tubeCurves2DImageFilter.cxx
  )

add_library( TubeTKTubeSegmentation 
  ${TubeTK_Base_TubeSegmentation_CXX_Files}
  )

target_link_libraries( TubeTKTubeSegmentation 
  ${ITK_LIBRARIES}
  )

#install(FILES 
  #${TubeTK_Base_TubeSegmentation_H_Files}
  #${TubeTK_Base_TubeSegmentation_H_Files}
  #DESTINATION include
  #)


