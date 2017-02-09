#add_executable("dune-PFASST" dune-PFASST.cc)
add_executable("FE_pfasst" FE_pfasst.cpp)
target_link_dune_default_libraries("FE_pfasst")
