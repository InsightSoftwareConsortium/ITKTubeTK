set(TubeTK_Base_Preprocessing_H_Files
  itkVesselEnhancingDiffusion2DImageFilter.h
)

set(TubeTK_Base_Preprocessing_TXX_Files
  itkVesselEnhancingDiffusion2DImageFilter.txx
)

install(FILES 
  ${TubeTK_Base_Preprocessing_H_Files}
  ${TubeTK_Base_Preprocessing_H_Files}
  DESTINATION include)