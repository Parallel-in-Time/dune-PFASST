


add_executable("FE_sdc_" FE_sdc.cpp)
target_link_dune_default_libraries("FE_sdc_")

add_executable("FE_mlsdc_" FE_mlsdc.cpp)
target_link_dune_default_libraries("FE_mlsdc_")

add_executable("FE_pfasst_" FE_pfasst.cpp)
target_link_dune_default_libraries("FE_pfasst_")
