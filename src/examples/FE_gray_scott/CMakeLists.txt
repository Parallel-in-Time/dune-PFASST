add_executable("FE_sdcGS" FE_sdcGS.cpp)
target_link_dune_default_libraries("FE_sdcGS")

add_executable("FE_mlsdcGS" FE_mlsdcGS.cpp)
target_link_dune_default_libraries("FE_mlsdcGS")

add_executable("FE_pfasstGS" FE_pfasstGS.cpp)
target_link_dune_default_libraries("FE_pfasstGS")
