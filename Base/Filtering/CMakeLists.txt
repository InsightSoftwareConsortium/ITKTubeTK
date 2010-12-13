set(TubeTK_Base_Filtering_H_Files
  itkAnisotropicCoherenceEnhancingDiffusionImageFilter.h
  itkAnisotropicDiffusionTensorFunction.h
  itkAnisotropicDiffusionTensorImageFilter.h
  itkAnisotropicEdgeEnhancementDiffusionImageFilter.h
  itkAnisotropicHybridDiffusionImageFilter.h
  itkStructureTensorRecursiveGaussianImageFilter.h
  itkSymmetricEigenVectorAnalysisImageFilter.h
  itkTubeBlurImageFunction.h
  itkTubeEnhancingDiffusion2DImageFilter.h
)

set(TubeTK_Base_Filtering_TXX_Files
  itkTubeEnhancingDiffusion2DImageFilter.txx
  itkStructureTensorRecursiveGaussianImageFilter.txx
  itkAnisotropicHybridDiffusionImageFilter.txx
  itkAnisotropicEdgeEnhancementDiffusionImageFilter.txx
  itkAnisotropicDiffusionTensorImageFilter.txx
  itkAnisotropicCoherenceEnhancingDiffusionImageFilter.txx
  itkAnisotropicDiffusionTensorFunction.txx
  itkTubeBlurImageFunction.txx
  itkTubeEnhancingDiffusion2DImageFilter.txx
)

IF( BUILD_TESTING )
  add_subdirectory( Testing )
ENDIF( BUILD_TESTING )

#install(FILES
  #${TubeTK_Base_Filtering_H_Files}
  #${TubeTK_Base_Filtering_H_Files}
  #DESTINATION include)
