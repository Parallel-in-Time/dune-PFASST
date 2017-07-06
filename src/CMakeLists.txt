add_subdirectory(examples)


if(dune-tnnmg_FOUND)
  add_executable("pfasst-tnnmg" dune-pfasst-tnnmg.cc)
  target_link_dune_default_libraries("pfasst-tnnmg")
  add_executable("pfasst-tnnmg-dg" dune-pfasst-tnnmg-dg.cc)
  target_link_dune_default_libraries("pfasst-tnnmg-dg")
endif()

