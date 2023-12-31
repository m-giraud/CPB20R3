dune_symlink_to_source_files(FILES "grids" "params.input")

# test gravity using different values for xi
add_executable(test_md_facet_1p1p_gravity_tpfa EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_md_facet_1p1p_gravity_tpfa
                           PUBLIC BULKTYPETAG=OnePBulkTpfa
                                  LOWDIMTYPETAG=OnePLowDimTpfa
                                  LOWDIMGRIDTYPE=Dune::FoamGrid<1,2>
                                  BULKGRIDTYPE=Dune::ALUGrid<2,2,Dune::cube,Dune::nonconforming>)

dumux_add_test(NAME test_md_facet_1p1p_gravity_xi1_tpfa
              LABELS multidomain
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TARGET test_md_facet_1p1p_gravity_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi1_bulk-00000.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi1_bulk-00001.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi1_lowdim-00000.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi1_lowdim-00001.vtp
                       --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_tpfa params.input")

dumux_add_test(NAME test_facet_1p1p_gravity_xi066_tpfa
              LABELS multidomain
              CMAKE_GUARD "( dune-foamgrid_FOUND AND dune-alugrid_FOUND )"
              TARGET test_md_facet_1p1p_gravity_tpfa
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS  --script fuzzy
                        --files ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi066_bulk-00000.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi066_bulk-00001.vtu
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi066_lowdim-00000.vtp
                                ${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_xi066_lowdim-00001.vtp
                        --command "${CMAKE_CURRENT_BINARY_DIR}/test_md_facet_1p1p_gravity_tpfa params.input  \
                                                               -FacetCoupling.Xi 0.66 -Vtk.OutputName test_md_facet_1p1p_gravity_xi066")
