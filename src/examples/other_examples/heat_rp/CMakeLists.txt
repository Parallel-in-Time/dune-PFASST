


add_executable("rsdc" FE_sdc.cpp)
target_link_dune_default_libraries("rsdc")

add_executable("rmlsdc_" FE_mlsdc.cpp)
target_link_dune_default_libraries("rmlsdc_")

add_executable("rpfasst_" FE_pfasst.cpp)
target_link_dune_default_libraries("rpfasst_")
