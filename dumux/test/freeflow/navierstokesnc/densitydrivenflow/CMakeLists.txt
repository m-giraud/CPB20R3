
dumux_add_test(NAME test_ff_stokes2c_densitydrivenflow
              SOURCES main.cc
              LABELS freeflow
              CMAKE_GUARD HAVE_UMFPACK
              COMMAND ${CMAKE_SOURCE_DIR}/bin/testing/runtest.py
              CMD_ARGS       --script fuzzy
                             --files ${CMAKE_SOURCE_DIR}/test/references/test_ff_stokes2c_densitydrivenflow-reference.vtu
                                     ${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_densitydrivenflow-00021.vtu
                             --command "${CMAKE_CURRENT_BINARY_DIR}/test_ff_stokes2c_densitydrivenflow params.input
                             -Problem.Name test_ff_stokes2c_densitydrivenflow"
                             --zeroThreshold {"X^Air_liq":1e-12}
                             --zeroThreshold {"x^Air_liq":1e-12}
                             --zeroThreshold {"velocity_liq \(m/s\)":1e-10}
                             --zeroThreshold {"deltaRho":1e-7})

dune_symlink_to_source_files(FILES "params.input")
