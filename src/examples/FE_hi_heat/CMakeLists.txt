
#add_executable("FE_hi_pfasst" FE_pfasst.cpp)
#target_link_dune_default_libraries("FE_hi_pfasst")

add_executable("FE_hi_sdc" FE_sdc.cpp)
target_link_dune_default_libraries("FE_hi_sdc")

add_executable("FE_hi_mlsdc" FE_mlsdc.cpp)
target_link_dune_default_libraries("FE_hi_mlsdc")
