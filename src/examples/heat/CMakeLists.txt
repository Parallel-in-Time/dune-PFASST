#add_executable("dune-PFASST" dune-PFASST.cc)
add_executable("dune-PFASST" vanilla_sdc.cpp)
target_link_dune_default_libraries("dune-PFASST")
