#add_executable("dune-PFASST" dune-PFASST.cc)
add_executable("pfasst" FE_pfasst.cpp)
target_link_dune_default_libraries("pfasst")

add_executable("sdc" FE_sdc.cpp)
target_link_dune_default_libraries("sdc")

add_executable("mlsdc" FE_mlsdc.cpp)
target_link_dune_default_libraries("mlsdc")
