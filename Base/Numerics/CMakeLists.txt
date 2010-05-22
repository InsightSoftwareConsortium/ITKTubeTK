set(TubeTK_Base_Numerics_H_Files
  itkNJetImageFunction.h
)

set(TubeTK_Base_Numerics_TXX_Files
  itkNJetImageFunction.txx
)

IF( BUILD_TESTING )
  add_subdirectory( Testing )
ENDIF( BUILD_TESTING )

#install(FILES 
  #${TubeTK_Base_Preprocessing_H_Files}
  #${TubeTK_Base_Preprocessing_H_Files}
  #DESTINATION include)

