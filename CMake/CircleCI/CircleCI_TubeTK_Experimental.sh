#!/bin/sh

MachineName=CircleCI_GitHub
BuildType=Release
CTestCommand=ctest
DashboardDir=/usr/src/

cd ${DashboardDir}
${CTestCommand} -D Experimental -D SITE_CTEST_MODE:STRING=Experimental -D SITE_BUILD_TYPE:STRING=${BuildType} -S /usr/src/ITKTubeTK/CMake/CircleCI/CircleCI_$1_Docker.cmake -V -VV -O CircleCI_$1_Docker.log
