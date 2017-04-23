


add_executable("FE_sdc" FE_sdc.cpp)
target_link_dune_default_libraries("FE_sdc")

add_executable("FE_mlsdc" FE_mlsdc.cpp)
target_link_dune_default_libraries("FE_mlsdc")
