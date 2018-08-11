


add_executable("istl_sdc" FE_sdc.cpp)
target_link_dune_default_libraries("istl_sdc")

add_executable("istl_mlsdc" FE_mlsdc.cpp)
target_link_dune_default_libraries("istl_mlsdc")

add_executable("istl_pfasst" FE_pfasst.cpp)
target_link_dune_default_libraries("istl_pfasst")
