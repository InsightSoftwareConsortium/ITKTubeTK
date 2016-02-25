
#
# TubeTK Maintainers:
#
# TubeTK base components required to build ITKTubeTK should be added to ITKTubeTK_DEPENDS list.
#
# * Components should be listed in topological order.
#
# * Name should be specified without the 'TubeTK' prefix.
#
# * All components should be listed: This will ensure all required folder will be added
# when building TubeTK as an ITK remote module.
#

set( ITKTubeTK_DEPENDS
  Common
  Numerics
  Filtering )

set( ITKTubeTK_LIBRARIES )
foreach( component ${ITKTubeTK_DEPENDS} )
  list( APPEND ITKTubeTK_LIBRARIES TubeTK${component} )
endforeach()
