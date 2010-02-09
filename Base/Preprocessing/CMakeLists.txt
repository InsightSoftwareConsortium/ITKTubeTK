set(TubeTK_Base_Preprocessing_H_Files
  itkTubeEnhancingDiffusion2DImageFilter.h
)

set(TubeTK_Base_Preprocessing_TXX_Files
  itkTubeEnhancingDiffusion2DImageFilter.txx
)

IF( BUILD_TESTING )
  add_subdirectory( Testing )
ENDIF( BUILD_TESTING )

#install(FILES 
  #${TubeTK_Base_Preprocessing_H_Files}
  #${TubeTK_Base_Preprocessing_H_Files}
  #DESTINATION include)
