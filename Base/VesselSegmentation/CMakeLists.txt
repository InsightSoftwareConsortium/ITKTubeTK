set(TubeTK_Base_VesselSegmentation_H_Files
  tubeCurves2DImageFilter.h
)

set(TubeTK_Base_VesselSegmentation_TXX_Files
)

set(TubeTK_Base_VesselSegmentation_CXX_Files
  tubeCurves2DImageFilter.cxx
)

add_library( TubeTKVesselSegmentation 
  ${TubeTK_Base_VesselSegmentation_CXX_Files}
)

target_link_libraries( TubeTKVesselSegmentation 
  ${ITK_LIBRARIES}
)

#install(FILES 
  #${TubeTK_Base_VesselSegmentation_H_Files}
  #${TubeTK_Base_VesselSegmentation_H_Files}
  #DESTINATION include)

#install(TARGETS TubeTKVesselSegmentation
  #DESTINATION lib)
