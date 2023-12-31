dune_symlink_to_source_files(FILES "params.input" "params_navierstokes.input" "params_convection.input" "params_conduction.input")

add_executable(test_ff_channel EXCLUDE_FROM_ALL main.cc)

dumux_add_test(NAME test_ff_stokes_channel
              TARGET test_ff_channel
              LABELS freeflow
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes_channel-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes_channel-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_channel params.input
                             -Problem.Name test_ff_stokes_channel")

dumux_add_test(NAME test_ff_navierstokes_channel
              TARGET test_ff_channel
              LABELS freeflow
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_navierstokes_channel-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_navierstokes_channel-00002.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_channel params_navierstokes.input
                             -Vtk.WriteFaceData 1
                             -Problem.Name test_ff_navierstokes_channel")

add_executable(test_ff_stokesni_channel EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_ff_stokesni_channel PUBLIC "NONISOTHERMAL=1")

dumux_add_test(NAME test_ff_stokesni_channel_convection
              TARGET test_ff_stokesni_channel
              LABELS freeflow
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokesni_channel_convection-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokesni_channel_convection-00006.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokesni_channel params_convection.input
                             -Problem.Name test_ff_stokesni_channel_convection")

dumux_add_test(NAME test_ff_stokesni_channel_conduction
              TARGET test_ff_stokesni_channel
              LABELS freeflow
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokesni_channel_conduction-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokesni_channel_conduction-00004.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokesni_channel params_conduction.input
                             -Problem.Name test_ff_stokesni_channel_conduction"
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-20})
