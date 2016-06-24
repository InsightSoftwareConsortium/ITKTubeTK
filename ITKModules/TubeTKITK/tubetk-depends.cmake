
#
# TubeTK Maintainers:
#
# TubeTK base components required to build TubeTKITK should be added to
#   TubeTKITK_DEPENDS list.
#
# * Components should be listed in topological order.
#
# * Name should be specified without the 'TubeTK' prefix.
#
# * All components should be listed: This will ensure all required folder
#     will be added when building TubeTK as an ITK remote module.
#

set( TubeTKITK_DEPENDS
  Common
  Numerics
  Filtering
  Segmentation  )

set( TubeTKITK_LIBRARIES )
foreach( component ${TubeTKITK_DEPENDS} )
  list( APPEND TubeTKITK_LIBRARIES TubeTK${component} )
endforeach()
