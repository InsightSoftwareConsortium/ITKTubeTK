#include "itkConfigure.h"
// Workaround for ITKTubeTK Python package builds
#if !defined(ITK_BUILD_SHARED_LIBS) && !defined(MinimalPathExtraction_EXPORT)
#  define MinimalPathExtraction_EXPORT
#  define MinimalPathExtraction_HIDDEN
#else
#include "MinimalPathExtractionExport.h"
#endif
