set(TESTS
  test-alugrid
  test-adaptivity
  test-fromtogridfactory
  test-structuredgridfactory
  test-twists
  test-hierarchic
  test-ball-projection
  test-backup-restore
)

foreach(TEST ${TESTS})
  dune_add_test(NAME ${TEST} SOURCES ${TEST}.cc LINK_LIBRARIES dunealugrid)
endforeach()

add_subdirectory(dgf)
