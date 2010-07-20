
include( ${CTEST_SCRIPT_DIRECTORY}/../../tubetk_config.cmake )

set( RUN_DASHBOARD_MODEL "Continuous" )

ctest_empty_binary_directory( "${SITE_BINARY_DIR}" )

###########################################################################
# run some "inside-the-loop" continuous scripts for a while
#
while(${CTEST_ELAPSED_TIME} LESS 56000)

  set(START_TIME ${CTEST_ELAPSED_TIME})

  ctest_run_script( 
    "${SITE_SCRIPT_DIR}/tubetk_build_test.cmake" )

  if( $ENV{TUBETK_CONTINUOUS_UPDATE} == 1 )

    if( SITE_CONTINUOUS_STYLECHECK )
      ctest_run_script( 
        "${SITE_SCRIPT_DIR}/tubetk_style.cmake" )
    endif( SITE_CONTINUOUS_STYLECHECK )
  
    if( SITE_CONTINUOUS_COVERAGE )
      ctest_run_script( 
        "${SCRIPT_DIR}/tubetk_coverage.cmake" )
    endif( SITE_CONTINUOUS_COVERAGE )
  
    if( SITE_CONTINUOUS_MEMORYCHECK )
      ctest_run_script( 
        "${SCRIPT_DIR}/tubetk_memory.cmake" )
    endif( SITE_CONTINUOUS_MEMORYCHECK )

  endif( $ENV{TUBETK_CONTINUOUS_UPDATE} == 1 )

  # loop no faster than once every 5 minutes
  ctest_sleep(${START_TIME} 300 ${CTEST_ELAPSED_TIME})

endwhile(${CTEST_ELAPSED_TIME} LESS 56000)

# do not run this script as a dashboard
set(CTEST_RUN_CURRENT_SCRIPT 0)

