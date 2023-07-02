#install headers
install(FILES
  computeintersection.hh
  computeintersection.cc
  conformingmerge.hh
  contactmerge.cc
  contactmerge.hh
  intersectionlist.hh
  merger.hh
  overlappingmerge.cc
  overlappingmerge.hh
  simplexintersection.cc
  standardmerge.hh
  DESTINATION include/dune/grid-glue/merging)

dune_add_library(dunegridglue
  standardmerge.cc
  ADD_LIBS ${DUNE_LIBS})
