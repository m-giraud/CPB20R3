add_subdirectory(staggered)

install(FILES
fluxvariables.hh
fluxvariablescache.hh
indices.hh
iofields.hh
localresidual.hh
model.hh
problem.hh
volumevariables.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/freeflow/navierstokes)
