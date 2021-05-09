find_package(VTK QUIET
 HINTS ./VTK-build ../VTK-build ../../VTK-build)
# HINTS are provided so that VTK is discovered in
#   build scripts used in Github Workflow actions.
#   See .github/workflows/build-test-package.yml
