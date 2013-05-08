#---------------------------------------------------------------
set(MODULE_NAME InteractiveConnectedComponentsUsingParzenPDFs)

#---------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
   ${MODULE_NAME}.py
   )

#---------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  )
