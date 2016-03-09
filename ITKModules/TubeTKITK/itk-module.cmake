set(DOCUMENTATION
  "ITK conformant classes available as an ITK external module")

itk_module( TubeTKITK
  DEPENDS
    ITKCommon
  EXCLUDE_FROM_DEFAULT
  DESCRIPTION
    "${DOCUMENTATION}"
)
