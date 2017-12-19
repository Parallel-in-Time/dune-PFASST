


add_executable("pde_sdc" FE_sdc.cpp)
target_link_dune_default_libraries("pde_sdc")

#add_executable("FE_mlsdc" FE_mlsdc.cpp)
#target_link_dune_default_libraries("FE_mlsdc")

#add_executable("FE_pfasst" FE_pfasst.cpp)
#target_link_dune_default_libraries("FE_pfasst")
