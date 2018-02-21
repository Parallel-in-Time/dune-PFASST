add_executable("FE_sdcGS_" FE_sdcGS.cpp)
target_link_dune_default_libraries("FE_sdcGS_")

add_executable("FE_mlsdcGS_" FE_mlsdcGS.cpp)
target_link_dune_default_libraries("FE_sdcGS_")

add_executable("FE_pfasst_" FE_pfasstGS.cpp)
target_link_dune_default_libraries("FE_pfasst_")

